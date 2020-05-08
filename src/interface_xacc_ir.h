#pragma once

#define DO_PRAGMA(x) _Pragma(#x)

#ifndef TODO
#define TODO(x) DO_PRAGMA(message("\033[30;43mTODO\033[0m - " #x))
#endif

#define XACC_QUAC_API __attribute__ ((visibility ("default"))) extern

// Interface to XACC IR: 
// Simulation mode: Circuit (gates) or Pulse (Hamiltonian)
typedef enum {
    CIRCUIT  = 0,
    PULSE = 1
} sim_mode;

typedef enum {
    NONE  = 0,
    MINIMAL = 1, // Important logging
    DEBUG = 2, // More logging
    DEBUG_DIAG = 3 // Very verbose
} log_verbosity;

typedef struct ComplexCoefficient {
    double real;
    double imag;
} ComplexCoefficient;

// Time-stepping data
typedef struct TSData {
    double time;
    int nbChannels;
    double* channelData;
    int nbPops;
    double* populations;
    double* pauliExpectations;
} TSData;

typedef struct XaccPulseChannelProvider PulseChannelProvider;



// ============================= Configuration API's =================================
// Pulse simulation initialization:
// Note: we *solve* the master equation using QuaC, not via Monte-Carlo method.
// Hence, we don't need to have the *shots* params.
XACC_QUAC_API int XACC_QuaC_Initialize(int in_nbQubit, const int* in_qbitDims);

// Clean-up any allocated resources
XACC_QUAC_API void XACC_QuaC_Finalize();

// Control logging verbosity
XACC_QUAC_API void XACC_QuaC_SetLogVerbosity(log_verbosity in_verboseConfig);

// Disable adaptive time-stepping (e.g. when time-dependent functions are oscillating)
XACC_QUAC_API void XACC_QuaC_DisableAdaptiveTimestepping();

// Set initial qubit population
XACC_QUAC_API void XACC_QuaC_SetInitialPopulation(int in_qubitIdx, double in_initialPopulation);

// Add a decay term (Lindblad)
XACC_QUAC_API void XACC_QuaC_AddQubitDecay(int in_qubitIdx, double in_kappa);

// Enable time-stepping data collection
// Specify the dt (e.g. sample duration) to save data
// The TS solver often has a much smaller dt which we don't want to get all data at each step.
XACC_QUAC_API void XACC_QuaC_EnableTimeSteppingMonitor(double in_dt);
// =========================================================================


// ============================= Hamiltonian construction API's =================================

// Adding a single-operator term to the Hamiltonian:
// (1) Time-independent term: 
// Syntax: coeff * ['X', 'Y', 'Z', 'I', 'SP', 'SM', 'O', 'N']_i
// 'SP' and 'SM' are the sigma plus and sigma minus operators.
// Coefficient is a complex parameter and this term can only act on 1 qubit.
XACC_QUAC_API void XACC_QuaC_AddConstHamiltonianTerm1(const char* in_op, int in_qubitIdx, ComplexCoefficient in_coeff);

// (2) Time-dependent term:
// Similar to (1) but has a time-dependent drive function (double -> double)
// Note: drive signal must have been *mixed* with LO, i.e. it is Re[d(t) * exp(-i * w_LO * t)] = d(t) * cos(w_LO * t)
XACC_QUAC_API void XACC_QuaC_AddTimeDependentHamiltonianTerm1(const char* in_op, int in_qubitIdx, int in_channelId, double in_coefficient);


// Adding a two-operator term to the Hamiltonian:
// (1) Time-independent term: 
XACC_QUAC_API void XACC_QuaC_AddConstHamiltonianTerm2(const char* in_op1, int in_qubitIdx1, const char* in_op2, int in_qubitIdx2, ComplexCoefficient in_coeff);

// (2) Time-dependent term:
XACC_QUAC_API void XACC_QuaC_AddTimeDependentHamiltonianTerm2(const char* in_op1, int in_qubitIdx1, const char* in_op2, int in_qubitIdx2, int in_channelId, double in_coefficient);

// ==============================================================================================


// ============================= Digital Gate API's =================================

// Add a U3 digital instruction to QuaC at a specific time 
// return 0 Success - operatation is added; 
// otherwise (>0) Failure - operatation cannot be added
XACC_QUAC_API int XACC_QuaC_AddDigitalInstructionU3(int in_qubitIdx, double in_theta, double in_phi, double in_lambda, double in_startTime);

// Add a digital CNOT gate:
// Note: we only support U3 and CNOT gates to be simulated in digital mode.
// (in case pulse library don't have the equivalent)
XACC_QUAC_API int XACC_QuaC_AddCnot(int in_ctrlIdx, int in_targetIdx, double in_startTime);
// ==================================================================================

// ============================= Simulation and Result API's =================================
 
// Run the Pulse simulation and return the expectation values:
// Returns the size of the result array. Caller needs to clean up. 
XACC_QUAC_API int XACC_QuaC_RunPulseSim(PulseChannelProvider* in_pulseDataProvider, double in_dt, double in_stopTime, int in_stepMax, double** out_result, int* out_nbSteps, TSData** out_timeSteppingData);

// Some methods to compute useful properties from the density matrix, e.g. to append as extra info to the XACC accelerator buffer
// Calculate the concurrence between 2 qubits
XACC_QUAC_API double XACC_QuaC_CalcConcurrence(int in_qubitIdx1, int in_qubitIdx2);

// Get a specific element from the density matrix
XACC_QUAC_API ComplexCoefficient XACC_QuaC_GetDensityMatrixElement(int in_row, int in_column);

// Get diagonal elements of the density matrix
XACC_QUAC_API int XACC_QuaC_GetDensityMatrixDiagElements(ComplexCoefficient** out_result);

// Calc. fidelity of the result density matrix against a reference 
XACC_QUAC_API double XACC_QuaC_CalcDensityMatrixFidelity(int in_size, ComplexCoefficient* in_refDm);
// ==============================================================================================
