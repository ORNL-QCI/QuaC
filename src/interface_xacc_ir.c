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
int nbQubits; 

bool g_wasInitialized = false;
int g_nbTimeDepChannels = 0;
bool g_enableTimeSteppingDataCollection = false;

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
// Register 8 channels
// Technically, we can register as many channels as we want here 
// (need to automate the macro expansion util to handle arbitrary number).
// These are just placeholder functions for calling into the Pulse controller.
typedef double channelFunctionType(double time);
REGISTER_N_DRIVE_CHANNELS(g_PulseDataProvider, 8);
// ***Temporary code for testing***
// This should be combined with the above REGISTER_N_DRIVE_CHANNELS macro.
channelFunctionType *g_channelFnArray[8] = { _DriveChannel0, _DriveChannel1, _DriveChannel2, _DriveChannel3, _DriveChannel4, _DriveChannel5, _DriveChannel6, _DriveChannel7 };

// Time-stepping monitor function
PetscErrorCode g_tsDefaultMonitorFunc(TS, PetscInt, PetscReal, Vec, void*);

int XACC_QuaC_Initialize(int in_nbQubit)
{
    // TODO: support passing params (e.g. qubit decay)
    QuaC_initialize(0, NULL);
    nbQubits = in_nbQubit;
    qubits  = malloc(in_nbQubit * sizeof(struct operator));    
    
    for (int i = 0; i < nbQubits; i++)
    {
        create_op(2, &qubits[i]);
    }
    
    // Allocate resources
    create_full_dm(&psi);
    // Initital state
    add_value_to_dm(psi, 0, 0, 1.0);

    // Create circuit
    create_circuit(&g_circuit, -1);
    return 0;
}

int XACC_QuaC_AddDigitalInstructionU3(int in_qubitIdx, double in_theta, double in_phi, double in_lambda, double in_startTime)
{
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

const char* XACC_QuaC_ExecuteCircuit(int in_argCount, char** in_args)
{
    start_circuit_at_time(&g_circuit, 0.0);
    time_step(psi, 0.0, g_timeMax, g_dt, g_stepsMax);
    // Debug
    int psi_dims = pow(2, nbQubits);
    print_psi(psi, psi_dims);
    return "";
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
    }
   
    _num_circuits = 0;
    _current_circuit = 0;
    g_gateStartTime = 0.0;
    QuaC_clear();
}

int XACC_QuaC_InitializePulseSim(int in_nbQubit, PulseChannelProvider* in_pulseDataProvider, const int* in_qbitDims)
{
    g_simulationMode = PULSE;
    g_PulseDataProvider = in_pulseDataProvider;
    if (!g_wasInitialized)
    {
        QuaC_initialize(0, NULL);
        g_wasInitialized = true;
    }
   
    nbQubits = in_nbQubit;
    qubits  = malloc(in_nbQubit * sizeof(struct operator));    
    
    for (int i = 0; i < nbQubits; i++)
    {
        int qubitDim = in_qbitDims[i];
        LOG_INFO("Qubit %d : Dim = %d \n", i, qubitDim);
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


void XACC_QuaC_AddConstHamiltonianTerm1(const char* in_op, int in_qubitIdx, ComplexCoefficient in_coeff)
{
    ASSERT_PULSE_MODE;
    ASSERT_QUBIT_INDEX(in_qubitIdx);
    
    LOG_INFO("H += (%lf + 1j*%lf)*%s%d\n", in_coeff.real, in_coeff.imag, in_op, in_qubitIdx);

    add_to_ham(in_coeff.real + in_coeff.imag * PETSC_i, GetQubitOperator(qubits[in_qubitIdx], in_op));
}


void XACC_QuaC_AddTimeDependentHamiltonianTerm1(const char* in_op, int in_qubitIdx, int in_channelId, double in_coefficient)
{
    ASSERT_PULSE_MODE;
    ASSERT_QUBIT_INDEX(in_qubitIdx);
    
    LOG_INFO("H += %lf * %s%d * Channel_%d (t)\n", in_coefficient, in_op, in_qubitIdx, in_channelId);

    add_to_ham_time_dep_with_coeff(in_coefficient, g_channelFnArray[in_channelId], 1, GetQubitOperator(qubits[in_qubitIdx], in_op));
    // Number of channels is the max of index + 1
    g_nbTimeDepChannels = (in_channelId + 1 > g_nbTimeDepChannels) ? in_channelId + 1 : g_nbTimeDepChannels;
}


void XACC_QuaC_AddConstHamiltonianTerm2(const char* in_op1, int in_qubitIdx1, const char* in_op2, int in_qubitIdx2, ComplexCoefficient in_coeff)
{
    ASSERT_PULSE_MODE;
    ASSERT_QUBIT_INDEX(in_qubitIdx1);
    ASSERT_QUBIT_INDEX(in_qubitIdx2);
    
    LOG_INFO("H += (%lf + 1j*%lf)*%s%d*%s%d\n", in_coeff.real, in_coeff.imag, in_op1, in_qubitIdx1, in_op2, in_qubitIdx2);

    add_to_ham_mult2(in_coeff.real + in_coeff.imag * PETSC_i, GetQubitOperator(qubits[in_qubitIdx1], in_op1), GetQubitOperator(qubits[in_qubitIdx2], in_op2));
}

void XACC_QuaC_AddTimeDependentHamiltonianTerm2(const char* in_op1, int in_qubitIdx1, const char* in_op2, int in_qubitIdx2, int in_channelId, double in_coefficient)
{
    ASSERT_PULSE_MODE;
    ASSERT_QUBIT_INDEX(in_qubitIdx1);
    ASSERT_QUBIT_INDEX(in_qubitIdx2);
    // Number of channels is the max of index + 1
    g_nbTimeDepChannels = (in_channelId + 1 > g_nbTimeDepChannels) ? in_channelId + 1 : g_nbTimeDepChannels;
    
    LOG_INFO("H += %lf * %s%d * %s%d * Channel_%d (t)\n", in_coefficient, in_op1, in_qubitIdx1, in_op2, in_qubitIdx2, in_channelId);
    add_to_ham_time_dep_with_coeff(in_coefficient, g_channelFnArray[in_channelId], 2, GetQubitOperator(qubits[in_qubitIdx1], in_op1), GetQubitOperator(qubits[in_qubitIdx2], in_op2));
}

int XACC_QuaC_RunPulseSim(double in_dt, double in_stopTime, int in_stepMax, double** out_result, int* out_nbSteps, TSData** out_timeSteppingData)
{
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

            g_timeSteppingData[g_nbStepCount].pauliExpectations[i*3] = PetscRealPart(expectX);
            g_timeSteppingData[g_nbStepCount].pauliExpectations[i*3 + 1] = PetscRealPart(expectY);
            g_timeSteppingData[g_nbStepCount].pauliExpectations[i*3 + 2] = PetscRealPart(expectZ);
        }
        else
        {
            g_timeSteppingData[g_nbStepCount].pauliExpectations[i*3] = 0.0;
            g_timeSteppingData[g_nbStepCount].pauliExpectations[i*3 + 1] = 0.0;
            g_timeSteppingData[g_nbStepCount].pauliExpectations[i*3 + 2] = 0.0;
        }     
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
