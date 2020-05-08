#include "interface_xacc_ir.h"
#include "quac.h"
#include <math.h>
#include "operators.h"
#include "solver.h"
#include "dm_utilities.h"
#include "quantum_gates.h"
#include "petsc.h"
#include <stdbool.h>
#include "macros.h"
#include "PulseControllerHandle.h"

// Global vars 
// Mode of simulation
// Default is circuit mode (quantum gates)
sim_mode g_simulationMode = CIRCUIT;

log_verbosity g_logVerboseLevel = MINIMAL;

void dispatchLog(log_verbosity in_level, const char* in_logFormat, ...) 
{
    // Skip all logs if NONE selected.
    if (g_logVerboseLevel == NONE)
    {
        return;
    }

    const bool shouldLog = (in_level <= g_logVerboseLevel);
    if (shouldLog)
    {
        va_list logArgs;
        va_start(logArgs, in_logFormat);
        vprintf(in_logFormat, logArgs);
        va_end(logArgs);      
    }
}

#define LOG_CRITICAL(format, ...) dispatchLog(MINIMAL, format, ##__VA_ARGS__)
#define LOG_INFO(format, ...) dispatchLog(DEBUG, format, ##__VA_ARGS__)
#define LOG_DEBUG(format, ...) dispatchLog(DEBUG_DIAG, format, ##__VA_ARGS__)

operator* qubits;
Vec psi;
circuit g_circuit;

int g_nbStepCount = 0;
TSData* g_timeSteppingData = NULL;

PulseChannelProvider* g_PulseDataProvider;

// Hard-coded for testing
PetscReal g_timeMax  = 9;
PetscReal g_dt = 0.01;
PetscInt g_stepsMax = 1000;

PetscReal gate_time_step = 1.0;
int nbQubits = 0; 
// Dimesition of the total density matrix:
// e.g. 2^nbQubits for qubit systems.
int densityMatrixDim = 1;

bool g_wasInitialized = false;
int g_nbTimeDepChannels = 0;
bool g_enableTimeSteppingDataCollection = false;
PetscReal g_monitorDt = 1.0;
PetscReal g_monitorTime = 0.0;


double g_gateStartTime = 0.0;

#define ASSERT_QUBIT_INDEX(qubitIdx) \
    if (qubitIdx > nbQubits) \
    { \
        printf("ERROR! Qubit index is out-of-range!\n");\
        exit(1);\
    }

#define ASSERT_PULSE_MODE \
    if (g_simulationMode != PULSE) \
    { \
        printf("ERROR! PULSE simulation mode has not been initialized!\n");\
        exit(1);\
    }

// Generate dummy channel drive functions (signature double(double)) which 
// just refer to g_PulseDataProvider to get the data.
// (each function declared here has an implied channel index)
// Technically, we can register as many channels as we want here 
// (need to automate the macro expansion util to handle arbitrary number).
// These are just placeholder functions for calling into the Pulse controller.
// Register 40 channels (~20 qubits, this should be the max that we can realistically handle) 
typedef double channelFunctionType(double time);
REGISTER_N_DRIVE_CHANNELS(g_PulseDataProvider, 40);

// Time-stepping monitor function
PetscErrorCode g_tsDefaultMonitorFunc(TS, PetscInt, PetscReal, Vec, void*);

int XACC_QuaC_AddDigitalInstructionU3(int in_qubitIdx, double in_theta, double in_phi, double in_lambda, double in_startTime)
{
    ASSERT_QUBIT_INDEX(in_qubitIdx);
    LOG_INFO("Add U3(%lf,%lf,%lf) q[%d] @ t = %lf \n", in_theta, in_phi, in_lambda, in_qubitIdx, in_startTime);
    circuit  circ;
    create_circuit(&circ, 1);
    add_gate_to_circuit(&circ, 0.0, U3, in_qubitIdx, in_theta, in_phi, in_lambda);
    if (in_startTime > g_gateStartTime)
    {
        g_gateStartTime = in_startTime;
    }
    else
    {
        g_gateStartTime += g_dt;
    }
    
    start_circuit_at_time(&circ, g_gateStartTime);
    return 0;
}

int XACC_QuaC_AddCnot(int in_ctrlIdx, int in_targetIdx, double in_startTime)
{
    ASSERT_QUBIT_INDEX(in_ctrlIdx);
    ASSERT_QUBIT_INDEX(in_targetIdx);
    
    LOG_INFO("Add CNOT(q[%d], q[%d]) @ t = %lf \n", in_ctrlIdx, in_targetIdx, in_startTime);
    circuit  circ;
    create_circuit(&circ, 1);
    add_gate_to_circuit(&circ, 0.0, CNOT, in_ctrlIdx, in_targetIdx); 
    if (in_startTime > g_gateStartTime)
    {
        g_gateStartTime = in_startTime;
    }
    else
    {
        g_gateStartTime += g_dt;
    }
    
    start_circuit_at_time(&circ, g_gateStartTime);
    return 0;
}

void XACC_QuaC_Finalize()
{
    for (int i=0; i< nbQubits; i++)
    {
        destroy_op(&qubits[i]);
    }
    
    free(qubits);   
    destroy_dm(psi);

    if (g_enableTimeSteppingDataCollection)
    {
        if (g_nbStepCount > 0)
        {
            // Free population data (at each time-step)
            for (int i = 0; i < g_nbStepCount; i++)
            {
                free(g_timeSteppingData[i].populations);
                free(g_timeSteppingData[i].channelData);
                free(g_timeSteppingData[i].pauliExpectations);
            }
            // Free the array of timestepping data-structures itself.
            free(g_timeSteppingData);
        }
        g_nbStepCount = 0;
    }
   
    _num_circuits = 0;
    _current_circuit = 0;
    g_gateStartTime = 0.0;
    g_enableTimeSteppingDataCollection = false;
    g_monitorTime = 0.0;
    QuaC_clear();
}

int XACC_QuaC_Initialize(int in_nbQubit, const int* in_qbitDims)
{
    g_simulationMode = PULSE;
    if (!g_wasInitialized)
    {
        QuaC_initialize(0, NULL);
        g_wasInitialized = true;
    }
   
    nbQubits = in_nbQubit;
    qubits  = malloc(in_nbQubit * sizeof(struct operator));    
    densityMatrixDim = 1;
    
    for (int i = 0; i < nbQubits; i++)
    {
        int qubitDim = in_qbitDims[i];
        LOG_INFO("Qubit %d : Dim = %d \n", i, qubitDim);
        // Update total density matrix dimension
        densityMatrixDim = densityMatrixDim * qubitDim;
        create_op(qubitDim, &qubits[i]);
        set_initial_pop(qubits[i], 0);
    }
     
    return 0;
}

void XACC_QuaC_SetLogVerbosity(log_verbosity in_verboseConfig)
{
    LOG_CRITICAL("Set logging level to %d.\n", in_verboseConfig);
    g_logVerboseLevel = in_verboseConfig;
}

void XACC_QuaC_AddQubitDecay(int in_qubitIdx, double in_kappa)
{
    ASSERT_PULSE_MODE;
    ASSERT_QUBIT_INDEX(in_qubitIdx);
    add_lin(in_kappa, qubits[in_qubitIdx]);
}

bool IsQubit(int in_idx)
{
    return qubits[in_idx]->my_levels == 2;
}

bool IsPauliOp(const char* in_op)
{
    return strcmp(in_op, "X") == 0 ||  strcmp(in_op, "Y") == 0 || strcmp(in_op, "Z") == 0;
}

operator GetQubitOperator(operator qubitOp, const char* in_op)
{
    if (strcmp(in_op, "SM") == 0) {
        // Sigma minus, i.e. itself
        return qubitOp;
    } else if (strcmp(in_op, "SP") == 0) {
        return (qubitOp)->dag;
    } else if (strcmp(in_op, "X") == 0) {
        return (qubitOp)->sig_x;
    } else if (strcmp(in_op, "Y") == 0) {
        return (qubitOp)->sig_y;
    } else if (strcmp(in_op, "Z") == 0) {
        return (qubitOp)->sig_z;
    } else if (strcmp(in_op, "I") == 0) {
        return (qubitOp)->eye;
    } else if (strcmp(in_op, "O") == 0 || strcmp(in_op, "N") == 0) {
        return (qubitOp)->n;
    } else {
        printf("ERROR! Unknown operator!\n");
        exit(1);
    }
}

// Pauli X, Y, Z decomposition for dim > 2 systems
struct PauliDecomp 
{
    PetscScalar op1_coeff;
    operator op1;
    PetscScalar op2_coeff;
    operator op2;
};

struct PauliDecomp getPauliDecomposition(const char* in_op, int in_qubitIdx)
{
    if (strcmp(in_op, "X") == 0) 
    {
        struct PauliDecomp result = { 1.0 , qubits[in_qubitIdx]->dag, 1.0, qubits[in_qubitIdx] };
        return result;
    }
    else if (strcmp(in_op, "Y") == 0) 
    {
        struct PauliDecomp result = { 1.0 * PETSC_i, qubits[in_qubitIdx]->dag, -1.0 * PETSC_i, qubits[in_qubitIdx] };
        return result;
    }
    else if (strcmp(in_op, "Z") == 0)
    {
        struct PauliDecomp result = { 1.0, qubits[in_qubitIdx]->eye, -2.0, qubits[in_qubitIdx]->n };
        return result;
    }
    printf("ERROR! Unknown operator!\n");
    exit(1);
}

void XACC_QuaC_AddConstHamiltonianTerm1(const char* in_op, int in_qubitIdx, ComplexCoefficient in_coeff)
{
    ASSERT_PULSE_MODE;
    ASSERT_QUBIT_INDEX(in_qubitIdx);
    
    LOG_INFO("H += (%lf + 1j*%lf)*%s%d\n", in_coeff.real, in_coeff.imag, in_op, in_qubitIdx);
    // If the system is Qubit (dimension = 2) or the operator is not a Pauli X, Y, Z
    // just process as normal
    if (IsQubit(in_qubitIdx) || !IsPauliOp(in_op))
    {
        add_to_ham(in_coeff.real + in_coeff.imag * PETSC_i, GetQubitOperator(qubits[in_qubitIdx], in_op));
    }
    else
    {
        struct PauliDecomp pairOps = getPauliDecomposition(in_op, in_qubitIdx);
        // Add the two terms associated with the Pauli ops
        add_to_ham((in_coeff.real + in_coeff.imag * PETSC_i) * pairOps.op1_coeff, pairOps.op1);
        add_to_ham((in_coeff.real + in_coeff.imag * PETSC_i) * pairOps.op2_coeff, pairOps.op2);
    }
}


void XACC_QuaC_AddTimeDependentHamiltonianTerm1(const char* in_op, int in_qubitIdx, int in_channelId, double in_coefficient)
{
    ASSERT_PULSE_MODE;
    ASSERT_QUBIT_INDEX(in_qubitIdx);
    
    LOG_INFO("H += %lf * %s%d * Channel_%d (t)\n", in_coefficient, in_op, in_qubitIdx, in_channelId);
    if (IsQubit(in_qubitIdx) || !IsPauliOp(in_op))
    {
        add_to_ham_time_dep_with_coeff(in_coefficient, g_channelFnArray[in_channelId], 1, GetQubitOperator(qubits[in_qubitIdx], in_op));
    }
    else
    {
        struct PauliDecomp pairOps = getPauliDecomposition(in_op, in_qubitIdx);
        // Add the two terms associated with the Pauli ops
        // TODO: support complex coefficients (Y operator) here
        add_to_ham_time_dep_with_coeff(in_coefficient * PetscRealPart(pairOps.op1_coeff), g_channelFnArray[in_channelId], 1, pairOps.op1);
        add_to_ham_time_dep_with_coeff(in_coefficient * PetscRealPart(pairOps.op2_coeff), g_channelFnArray[in_channelId], 1, pairOps.op2);
    }

    // Number of channels is the max of index + 1
    g_nbTimeDepChannels = (in_channelId + 1 > g_nbTimeDepChannels) ? in_channelId + 1 : g_nbTimeDepChannels;
}


void XACC_QuaC_AddConstHamiltonianTerm2(const char* in_op1, int in_qubitIdx1, const char* in_op2, int in_qubitIdx2, ComplexCoefficient in_coeff)
{
    ASSERT_PULSE_MODE;
    ASSERT_QUBIT_INDEX(in_qubitIdx1);
    ASSERT_QUBIT_INDEX(in_qubitIdx2);
    
    LOG_INFO("H += (%lf + 1j*%lf)*%s%d*%s%d\n", in_coeff.real, in_coeff.imag, in_op1, in_qubitIdx1, in_op2, in_qubitIdx2);
    
    // Handle case where it's a product of two operators on the same qubit sub-system
    // i.e. not a Kronecker product
    // e.g., N1*N1 etc.
    // !!IMPORTANT!! QuaC solver will not be able to handle this case
    // add_to_ham_mult2 assumes two operators act on different qubit subspaces.
    if (in_qubitIdx1 == in_qubitIdx2)    
    {
        // Handle a product term between two operators on the same qubit subspace
        Mat tempMat;
        // Construct the product matrix
        //combine_ops_to_mat(&tempMat, 2, GetQubitOperator(qubits[in_qubitIdx1], in_op1), GetQubitOperator(qubits[in_qubitIdx1], in_op2));
        mult_ops_as_mat(&tempMat, GetQubitOperator(qubits[in_qubitIdx1], in_op1), GetQubitOperator(qubits[in_qubitIdx1], in_op2));
        PetscInt rows, cols;
        MatGetSize(tempMat, &rows, &cols);
        LOG_INFO("Product operator matrix size = %ld, %ld\n", rows, cols);
        
        add_mat_to_ham(in_coeff.real + in_coeff.imag * PETSC_i, qubits[in_qubitIdx1], tempMat);
        MatDestroy(&tempMat);
        return;
    }    
    
    if ((IsQubit(in_qubitIdx1) || !IsPauliOp(in_op1)) && (IsQubit(in_qubitIdx2) || !IsPauliOp(in_op2)))
    {
        add_to_ham_mult2(in_coeff.real + in_coeff.imag * PETSC_i, GetQubitOperator(qubits[in_qubitIdx1], in_op1), GetQubitOperator(qubits[in_qubitIdx2], in_op2));
    }
    else
    {
        bool q1NeedTransform = !IsQubit(in_qubitIdx1) && IsPauliOp(in_op1); 
        bool q2NeedTransform = !IsQubit(in_qubitIdx2) && IsPauliOp(in_op2); 
        if (q1NeedTransform && !q2NeedTransform)
        {
            struct PauliDecomp pairOps = getPauliDecomposition(in_op1, in_qubitIdx1);
            add_to_ham_mult2((in_coeff.real + in_coeff.imag * PETSC_i) * pairOps.op1_coeff, pairOps.op1, GetQubitOperator(qubits[in_qubitIdx2], in_op2));
            add_to_ham_mult2((in_coeff.real + in_coeff.imag * PETSC_i) * pairOps.op2_coeff, pairOps.op1, GetQubitOperator(qubits[in_qubitIdx2], in_op2));
        }
        else if (!q1NeedTransform && q2NeedTransform)
        {
            struct PauliDecomp pairOps = getPauliDecomposition(in_op2, in_qubitIdx2);
            add_to_ham_mult2((in_coeff.real + in_coeff.imag * PETSC_i) * pairOps.op1_coeff, GetQubitOperator(qubits[in_qubitIdx1], in_op1), pairOps.op1);
            add_to_ham_mult2((in_coeff.real + in_coeff.imag * PETSC_i) * pairOps.op2_coeff, GetQubitOperator(qubits[in_qubitIdx1], in_op1), pairOps.op1);
        }
        else
        {
            struct PauliDecomp pairOps1 = getPauliDecomposition(in_op1, in_qubitIdx1);
            struct PauliDecomp pairOps2 = getPauliDecomposition(in_op2, in_qubitIdx2);
            // Decompose both Pauli ops
            add_to_ham_mult2((in_coeff.real + in_coeff.imag * PETSC_i) * pairOps1.op1_coeff * pairOps2.op1_coeff, pairOps1.op1, pairOps2.op1);
            add_to_ham_mult2((in_coeff.real + in_coeff.imag * PETSC_i) * pairOps1.op1_coeff * pairOps2.op2_coeff, pairOps1.op1, pairOps2.op2);
            
            add_to_ham_mult2((in_coeff.real + in_coeff.imag * PETSC_i) * pairOps1.op2_coeff * pairOps2.op1_coeff, pairOps1.op2, pairOps2.op1);
            add_to_ham_mult2((in_coeff.real + in_coeff.imag * PETSC_i) * pairOps1.op2_coeff * pairOps2.op2_coeff, pairOps1.op2, pairOps2.op2);
        }
    }
}

void XACC_QuaC_AddTimeDependentHamiltonianTerm2(const char* in_op1, int in_qubitIdx1, const char* in_op2, int in_qubitIdx2, int in_channelId, double in_coefficient)
{
    ASSERT_PULSE_MODE;
    ASSERT_QUBIT_INDEX(in_qubitIdx1);
    ASSERT_QUBIT_INDEX(in_qubitIdx2);
    // Number of channels is the max of index + 1
    g_nbTimeDepChannels = (in_channelId + 1 > g_nbTimeDepChannels) ? in_channelId + 1 : g_nbTimeDepChannels;
    
    LOG_INFO("H += %lf * %s%d * %s%d * Channel_%d (t)\n", in_coefficient, in_op1, in_qubitIdx1, in_op2, in_qubitIdx2, in_channelId);
    if ((IsQubit(in_qubitIdx1) || !IsPauliOp(in_op1)) && (IsQubit(in_qubitIdx2) || !IsPauliOp(in_op2)))
    {
        add_to_ham_time_dep_with_coeff(in_coefficient, g_channelFnArray[in_channelId], 2, GetQubitOperator(qubits[in_qubitIdx1], in_op1), GetQubitOperator(qubits[in_qubitIdx2], in_op2));
    }
    else
    {
        bool q1NeedTransform = !IsQubit(in_qubitIdx1) && IsPauliOp(in_op1); 
        bool q2NeedTransform = !IsQubit(in_qubitIdx2) && IsPauliOp(in_op2); 
        if (q1NeedTransform && !q2NeedTransform)
        {
            struct PauliDecomp pairOps = getPauliDecomposition(in_op1, in_qubitIdx1);
            add_to_ham_time_dep_with_coeff(in_coefficient * PetscRealPart(pairOps.op1_coeff), g_channelFnArray[in_channelId], 2, pairOps.op1, GetQubitOperator(qubits[in_qubitIdx2], in_op2));
            add_to_ham_time_dep_with_coeff(in_coefficient * PetscRealPart(pairOps.op2_coeff), g_channelFnArray[in_channelId], 2, pairOps.op2, GetQubitOperator(qubits[in_qubitIdx2], in_op2));
        }
        else if (!q1NeedTransform && q2NeedTransform)
        {
            struct PauliDecomp pairOps = getPauliDecomposition(in_op2, in_qubitIdx2);
            add_to_ham_time_dep_with_coeff(in_coefficient * PetscRealPart(pairOps.op1_coeff), g_channelFnArray[in_channelId], 2, GetQubitOperator(qubits[in_qubitIdx1], in_op1), pairOps.op1);
            add_to_ham_time_dep_with_coeff(in_coefficient * PetscRealPart(pairOps.op1_coeff), g_channelFnArray[in_channelId], 2, GetQubitOperator(qubits[in_qubitIdx1], in_op1), pairOps.op2);
        }
        else
        {
            struct PauliDecomp pairOps1 = getPauliDecomposition(in_op1, in_qubitIdx1);
            struct PauliDecomp pairOps2 = getPauliDecomposition(in_op2, in_qubitIdx2);

            add_to_ham_time_dep_with_coeff(in_coefficient * PetscRealPart(pairOps1.op1_coeff) * PetscRealPart(pairOps2.op1_coeff), g_channelFnArray[in_channelId], 2, 
                pairOps1.op1, pairOps2.op1);
            
            add_to_ham_time_dep_with_coeff(in_coefficient * PetscRealPart(pairOps1.op1_coeff) * PetscRealPart(pairOps2.op2_coeff), g_channelFnArray[in_channelId], 2, 
                pairOps1.op1, pairOps2.op2);

            add_to_ham_time_dep_with_coeff(in_coefficient * PetscRealPart(pairOps1.op2_coeff) * PetscRealPart(pairOps2.op1_coeff), g_channelFnArray[in_channelId], 2, 
                pairOps1.op2, pairOps2.op1);
            add_to_ham_time_dep_with_coeff(in_coefficient * PetscRealPart(pairOps1.op2_coeff) * PetscRealPart(pairOps2.op2_coeff), g_channelFnArray[in_channelId], 2, 
                pairOps1.op2, pairOps2.op2);
        }
    }
}

int XACC_QuaC_RunPulseSim(PulseChannelProvider* in_pulseDataProvider, double in_dt, double in_stopTime, int in_stepMax, double** out_result, int* out_nbSteps, TSData** out_timeSteppingData)
{
    g_PulseDataProvider = in_pulseDataProvider;
    g_dt = in_dt;
    g_timeMax = in_stopTime;
    g_stepsMax = in_stepMax;
    
    // Allocate resources
    create_full_dm(&psi);
    set_dm_from_initial_pop(psi);

    if (g_enableTimeSteppingDataCollection)
    {
        set_ts_monitor(g_tsDefaultMonitorFunc);
    
        // Allocate an ample array for TS data.
        // Note: the real data (channels, populations, etc.) is allocated separarely, 
        // hence doesn't bloat the memory.
        g_timeSteppingData =  malloc(2 * g_stepsMax * sizeof(TSData));
    }

    time_step(psi, 0.0, g_timeMax, g_dt, g_stepsMax);
    // Returns the population for each qubit
    int nbResults = get_num_populations();
    *out_result = malloc(nbResults * sizeof(double));
    get_populations(psi, &(*out_result)); 
    
    if (g_enableTimeSteppingDataCollection)
    {    
        *out_nbSteps = g_nbStepCount;
        *out_timeSteppingData = g_timeSteppingData;
    }
    else
    {
        *out_nbSteps = 0;
    }
    

    return nbResults;
}

PetscErrorCode g_tsDefaultMonitorFunc(TS ts, PetscInt step, PetscReal time, Vec dm, void *ctx)
{
    if (time > g_monitorTime)
    {
        int num_pop = get_num_populations();
        double *populations;
        populations = malloc(num_pop*sizeof(double));
        get_populations(dm, &populations);
        PetscScalar expectX, expectY, expectZ;

        g_timeSteppingData[g_nbStepCount].time = time;
        g_timeSteppingData[g_nbStepCount].nbPops = num_pop;
        g_timeSteppingData[g_nbStepCount].populations = malloc(num_pop * sizeof(double));
        g_timeSteppingData[g_nbStepCount].nbChannels = g_nbTimeDepChannels;
        g_timeSteppingData[g_nbStepCount].channelData = malloc(g_nbTimeDepChannels * sizeof(double));
        g_timeSteppingData[g_nbStepCount].pauliExpectations = malloc(nbQubits * 3 * sizeof(double));
        memcpy(g_timeSteppingData[g_nbStepCount].populations, populations, num_pop * sizeof(double));

        for (int i = 0; i < g_nbTimeDepChannels; ++i)
        {
            g_timeSteppingData[g_nbStepCount].channelData[i] = g_channelFnArray[i](time);
        }

        for (int i = 0; i < nbQubits; ++i)
        {
            if (qubits[i]->my_levels == 2)
            {
                // Pauli-X
                get_expectation_value(dm, &expectX, 1, qubits[i]->sig_x);
                // Pauli-Y
                get_expectation_value(dm, &expectY, 1, qubits[i]->sig_y);
                // Pauli-Z
                get_expectation_value(dm, &expectZ, 1, qubits[i]->sig_z);
            }
            else
            {
                // Use the *custom* extended Pauli expectation calculation for non-qubit sub-systems.
                // Pauli-X
                expectX = get_expectation_value_Pauli_ext(dm, qubits[i]->sig_x);
                // Pauli-Y
                expectY = get_expectation_value_Pauli_ext(dm, qubits[i]->sig_y);
                // Pauli-Z
                expectZ = get_expectation_value_Pauli_ext(dm, qubits[i]->sig_z);
            }

            g_timeSteppingData[g_nbStepCount].pauliExpectations[i*3] = PetscRealPart(expectX);
            g_timeSteppingData[g_nbStepCount].pauliExpectations[i*3 + 1] = PetscRealPart(expectY);
            g_timeSteppingData[g_nbStepCount].pauliExpectations[i*3 + 2] = PetscRealPart(expectZ);     
        }

        if (nid==0)
        {
            // Debug:
            LOG_DEBUG(">> t = %e: ", time);
            for(int i = 0; i < num_pop; i++)
            {
                LOG_DEBUG("%e ", populations[i]);
            }
            LOG_DEBUG("\n");
        }

        free(populations);
        g_nbStepCount++;
        // Update the TS monitor time
        // We don't record all time steps of the solver
        g_monitorTime = g_monitorTime + g_monitorDt;
    }

    PetscFunctionReturn(0);
}

void XACC_QuaC_SetInitialPopulation(int in_qubitIdx, double in_initialPopulation)
{
    ASSERT_QUBIT_INDEX(in_qubitIdx);
    set_initial_pop(qubits[in_qubitIdx], in_initialPopulation);
}

void XACC_QuaC_DisableAdaptiveTimestepping()
{
    _disable_adaptive_ts = 1;
}

double XACC_QuaC_CalcConcurrence(int in_qubitIdx1, int in_qubitIdx2)
{
    ASSERT_QUBIT_INDEX(in_qubitIdx1);
    ASSERT_QUBIT_INDEX(in_qubitIdx2);
    // Partial trace: only keep the 2 requested qubits
    Vec tmpPsi;
    create_dm(&tmpPsi, 4);
    partial_trace_keep(psi, tmpPsi, 2, qubits[in_qubitIdx1], qubits[in_qubitIdx2]);
    
    // Calculate the bipartite concurrence
    double concurrenceResult = 0.0;
    get_bipartite_concurrence(tmpPsi, &concurrenceResult);

    destroy_dm(tmpPsi);

    return concurrenceResult;
}

ComplexCoefficient XACC_QuaC_GetDensityMatrixElement(int in_row, int in_column)
{
    PetscScalar dmElement;
    get_dm_element(psi, in_row, in_column, &dmElement);
    ComplexCoefficient result;
    result.real = PetscRealPart(dmElement);
    result.imag = PetscImaginaryPart(dmElement);
    return result;
}

int XACC_QuaC_GetDensityMatrixDiagElements(ComplexCoefficient** out_result)
{
    *out_result = malloc(densityMatrixDim * sizeof(ComplexCoefficient));
    // Get all diagonal elements of the density matrix
    for (int i = 0; i < densityMatrixDim; ++i)
    {
        (*out_result)[i] = XACC_QuaC_GetDensityMatrixElement(i, i);
    }
    
    return densityMatrixDim;
}

void XACC_QuaC_EnableTimeSteppingMonitor(double in_dt)
{
    g_enableTimeSteppingDataCollection = true;
    g_monitorDt = in_dt;
}

double XACC_QuaC_CalcDensityMatrixFidelity(int in_size, ComplexCoefficient* in_refDm)
{
    double resultFidelity = 0.0;
    // Make sure the dimension matched
    if (in_size == densityMatrixDim*densityMatrixDim)
    {
        Vec refDm;
        // Create a reference DM for fidelity calculations 
        create_dm(&refDm, densityMatrixDim);
        for (int rowIdx = 0; rowIdx < densityMatrixDim; ++rowIdx)
        {
            for (int colIdx = 0; colIdx < densityMatrixDim; ++colIdx)
            {
                int index = rowIdx*densityMatrixDim + colIdx;
                ComplexCoefficient val = in_refDm[index];
                add_value_to_dm(refDm, rowIdx, colIdx, val.real + val.imag * PETSC_i);
            }
        }
        assemble_dm(refDm);

        // Calculate the fidelity
        get_fidelity(psi, refDm, &resultFidelity);

        destroy_dm(refDm);
    }

    return resultFidelity;
}
