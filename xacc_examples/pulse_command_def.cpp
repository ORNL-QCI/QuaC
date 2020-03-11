#include "xacc.hpp"
#include <cassert>
#include "Pulse.hpp"
#include "xacc_service.hpp"


// ============================================================================
// This example demonstrates the backend configs loading, e.g. from IBM JSON file
// i.e., the Hamiltonian will be parsed directly from the Json file.
// ============================================================================
int main (int argc, char** argv) {
    // Initialize the XACC Framework
    xacc::Initialize(argc, argv);
    // IBM 20-qubit device Poughkeepsie backend JSON
    const std::string jsonConfigFile = std::string(GATEIR_TEST_FILE_DIR) + "/test_backends.json";
    
    // By providing a full backend config, the XACC-Quac accelerator will:
    // (1) Parse the system Hamiltonian from the Json
    // (2) Parse all the pulses that are defined (by samples) and all the so-called 'cmd-def''s
    // which are essentially sequence of pulses to realize a gate.
    auto quaC = xacc::getAccelerator("QuaC", 
        { std::make_pair("config-json-path", jsonConfigFile) });    

    // At this point, the whole pulse library has been imported.
    // For example, we can retrieve the pulse definition for a X(q[0]) gate for this device:
    auto x0 = xacc::ir::asComposite(xacc::getContributedService<xacc::Instruction>("pulse::x_0"));
    // For the Poughkeepsie device, X gate is just a single PI pulse named "Xp_d0_ddeb" (specific for qubit #0, i.e., d0 channel)
    assert(x0->nInstructions() == 1);
    assert(x0->getInstruction(0)->name() == "Xp_d0_ddeb");
    // Print out the definition
    std::cout << "X(q[0]): \n" << x0->toString() << "\n";
    
    // Similarly, CNOT (CX) gate is also defined in terms of pulses (and frame changes) 
    // which we have loaded from the JSON file. 
    // e.g. CNOT(q[0], q[1])
    auto cx01 = xacc::ir::asComposite(xacc::getContributedService<xacc::Instruction>("pulse::cx_0_1"));
    // For the Poughkeepsie device, CX gate is constructed from 10 commands (7 pulses and 3 frame changes)
    assert(cx01->nInstructions() == 10);
    
    // Print out the definition
    std::cout << "CNOT(q[0], q[1]): \n" << cx01->toString() << "\n";
    
    // Finalize the XACC Framework
    xacc::Finalize();

    return 0;
}