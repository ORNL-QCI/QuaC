
#include "xacc.hpp"
#include "Pulse.hpp"
#include "PulseSystemModel.hpp"
#include "PulseGen.hpp"

// Based on Eugene's Cross_Resonant_Interaction IPython Notebook
int main (int argc, char** argv) {
	// Initialize the XACC Framework
	xacc::Initialize(argc, argv);
    
    const auto runSim = [](bool in_controlVal) {
        auto systemModel = std::make_shared<QuaC::PulseSystemModel>();
        // Target Resonator's Reference Frame + RWA
        const std::string hamiltonianJson = R"#(
            {
                "description": "Two-qutrit Hamiltonian. Rotating frame @ target qubit frame",
                "h_latex": "",
                "h_str": ["(w_0-w_1)*O0", "d*O0*(O0-I0)", "d*O1*(O1-I1)", "J*(SP0*SM1 + SM0*SP1)", "O*(SM0 + SP0)||D0"],
                "osc": {},
                "qub": {
                    "0": 3,
                    "1": 3
                },
                "vars": {
                    "w_0": 5.114,
                    "w_1": 4.914,
                    "d": -0.33,
                    "J": 0.004,
                    "O": 0.060
                }
            }
        )#";

        if (!systemModel->loadHamiltonianJson(hamiltonianJson))
        {
            return;
        }
        
        BackendChannelConfigs channelConfigs;
        channelConfigs.dt = 5.0;
        // LO freq = 0.0 (rotating frame)
        channelConfigs.loFregs_dChannels.emplace_back(0.0);
        // Run to t = 5000 (equal the notebook)
        const size_t nbSamples = 1000;
        // Constant drive
        channelConfigs.addOrReplacePulse("square", QuaC::SquarePulse(nbSamples));        
        systemModel->setChannelConfigs(channelConfigs);    
        // Set the first qutrit (control) to |1> if needed 
        if (in_controlVal)
        {
            systemModel->setQubitInitialPopulation(0, 1.0);
        }
        // Run the simulation with time-stepping data collection
        auto quaC = xacc::getAccelerator("QuaC", { std::make_pair("system-model", systemModel), std::make_pair("logging-period", 0.1) });    
        auto qubitReg = xacc::qalloc(2);    

        auto provider = xacc::getIRProvider("quantum");
        auto compositeInst = provider->createComposite("test_pulse");
        auto pulseInstD0 = std::make_shared<xacc::quantum::Pulse>("square", "d0");

        compositeInst->addInstruction(pulseInstD0);
        // Run the Pulse simulation with the Hamiltonian provided
        quaC->execute(qubitReg, compositeInst);
    };
    
    // First run: control qutrit = 0 
    runSim(false);
    
    // Second run: control qutrit = 1 
    runSim(true);

    // Check the Console log to see the timestamped CSV files
    // that were exported.

    // Finalize the XACC Framework
	xacc::Finalize();

	return 0;
}