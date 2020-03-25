#include "xacc.hpp"
#include <cassert>
#include "Pulse.hpp"
#include "xacc_service.hpp"
#include "PulseSystemModel.hpp"

// ============================================================================
// Using cmd-def's of single-qubit gates from IBM Poughkeepsie backend 
// ============================================================================
int main (int argc, char** argv) {
	// Initialize the XACC Framework
	xacc::Initialize(argc, argv);
    
    // Connection Topology:
    // In this file, we try to find/calibrate missing Hamiltonian coefficients
    // such that the provided pulse library works reasonably well.
    // (e.g. the CNOT pulse cmd-def work as intended)
    // We targeted a subset of these qubits.
    // 0  --- 1 ---  2 ---  3 ---  4
    // |                           |
    // 5  --- 6 ---  7 ---  8 ---  9
    // |             |              |    
    // 10 --- 11 --- 12 --- 13 --- 14
    // |                           |
    // 15 --- 16 --- 17 --- 18 --- 19

    // For this one, we try to calibrate omegad0 and omegad1 (driving coefficient on D channels of each qubit)
    // and jq0q1 (coupling between Q0 and Q1).
    // Note that: we ignored any further couplings (1-2, 0-5) for now (considered as errors)
    const auto createHamiltonianJson = [](double in_omega0, double in_omega1, double in_j0j1) {
        const std::string hamiltonianJsonTmpl = R"(
        {
            "description": "Qubits are modelled as a two level system.\n",
            "h_str": ["_SUM[i,0,1,wq{i}/2*(I{i}-Z{i})]", "_SUM[i,0,1,omegad{i}*X{i}||D{i}]", "omegad1*X0||U0", "omegad5*X0||U1", "omegad0*X1||U2", "omegad2*X1||U3", "jq0q1*Sp0*Sm1", "jq0q1*Sm0*Sp1"],
            "osc": {},
            "qub": {
                "0": 2,
                "1": 2
            },
            "vars": {
                "omegad0": {{omegad0}},
                "omegad1": {{omegad1}}, 
                "omegad5": 0.0,
                "omegad2": 0.0,
                "wq0": 30.91270129264568,
                "wq1": 30.36010168900955,
                "jq0q1": {{jq0q1}}
            }
        })";

        const std::string omega0Placeholder = "{{omegad0}}";
        const std::string omega1Placeholder = "{{omegad1}}";
        const std::string j0j1Placeholder = "{{jq0q1}}";
        std::string hamiltonianJson = hamiltonianJsonTmpl;
        
        hamiltonianJson.replace(hamiltonianJson.find(omega0Placeholder), omega0Placeholder.length(), std::to_string(in_omega0));
        hamiltonianJson.replace(hamiltonianJson.find(omega1Placeholder), omega1Placeholder.length(), std::to_string(in_omega1));
        hamiltonianJson.replace(hamiltonianJson.find(j0j1Placeholder), j0j1Placeholder.length(), std::to_string(in_j0j1));
        
        return hamiltonianJson;
    };
    
    std::vector<double> loFreqs { 4.919909215047782, 4.83196025657847 };

    // U channels: cross-resonance frequencies.
    // Note: we disable 0-5 and 1-2 cross resonance (omegad5 and omegad2 are left at zero, we don't try to set these)
    // Just focus on the Q0-Q1 pair.
    std::vector<double> loFreqs_uChannels { 4.83196025657847, 4.957290417531115, 4.919909215047782, 4.940451402139736 };
    
    // The cmd-def for CNOT(Q[0], Q[1]) involves U10 channel, which is the cross-resonance channel on Q5.
    // Although we ignore Q5 in this calibration procedure, need to set U channel freqs 
    // (hence allocate channels in the channel controller)
    // Add extra 7 U channels (unused)
    for (int i = 0; i < 7; ++i)
    {
        loFreqs_uChannels.emplace_back(0.0);
    }

    // Run a composite instruction (lowered to pulse using the IBM pulse library) with a set of parameters
    // returns the diagonal element of the density matrix.
    // To be used as the kernel of an optimization procedure.
    const auto runSim = [&](const std::shared_ptr<xacc::CompositeInstruction>& in_program, double in_omega0, double in_omega1, double in_j0j1) {
        const std::string hamiltonianJson = createHamiltonianJson(in_omega0, in_omega1, in_j0j1);
        auto systemModel = std::make_shared<QuaC::PulseSystemModel>();
        const bool loadOk = systemModel->loadHamiltonianJson(hamiltonianJson);
        assert(loadOk);
        
        BackendChannelConfigs channelConfigs;
        channelConfigs.dt = 3.5555555555555554;
        channelConfigs.loFregs_dChannels = loFreqs;
        channelConfigs.loFregs_uChannels = loFreqs_uChannels;
        systemModel->setChannelConfigs(channelConfigs);    
        
        const int NB_SHOTS = 10000;
        auto quaC = xacc::getAccelerator("QuaC", { std::make_pair("system-model", systemModel), 
                                                    std::make_pair("shots", NB_SHOTS) });
        
        
        // Contribute all cmd-defs in the backend Json as XACC instruction.
        // This will activate Gates -> Pulses decomposition when simulating the circuit.
        static bool contributed = false;
        if (!contributed)
        {
            const std::string jsonConfigFile = std::string(GATEIR_TEST_FILE_DIR) + "/test_backends.json";
            std::ifstream backendFile(jsonConfigFile);
            std::string jjson((std::istreambuf_iterator<char>(backendFile)), std::istreambuf_iterator<char>());
            quaC->contributeInstructions(jjson);    
            contributed = true;
        }
        
        
        auto qubitReg = xacc::qalloc(2);
        
        // Run the Pulse simulation with the Hamiltonian provided
        quaC->execute(qubitReg, in_program);
        qubitReg->print();
        return qubitReg->getInformation("DensityMatrixDiags").as<std::vector<double>>();
    };

    // Create a programs using IR to calibrate
    auto provider = xacc::getIRProvider("quantum");
    auto x0 = provider->createInstruction("X", {0});
    auto x1 = provider->createInstruction("X", {1});
    auto h0 = provider->createInstruction("H", {0});
    auto h1 = provider->createInstruction("H", {1});
    auto cnot = provider->createInstruction("CNOT", {0, 1});
    
    auto compositeXX = provider->createComposite("test_XX");
    compositeXX->addInstructions({x0,x1});
    // Expected result (diagonal DM elements)
    std::vector<double> expectedCompositeXX { 0.0, 0.0, 0.0, 1.0 };

    auto compositeHH = provider->createComposite("test_HH");
    compositeHH->addInstructions({h0,h1});
    std::vector<double> expectedCompositeHH { 0.25, 0.25, 0.25, 0.25 };
    
    auto compositeXCNOT = provider->createComposite("test_XCNOT");
    compositeXCNOT->addInstructions({x0, cnot});
    std::vector<double> expectedCompositeXCNOT { 0.0, 0.0, 0.0, 1.0 };
    
    auto compositeHCNOT = provider->createComposite("test_HCNOT");
    compositeHCNOT->addInstructions({h0, cnot});
    std::vector<double> expectedCompositeHCNOT { 0.5, 0.0, 0.0, 0.5 };


    xacc::OptFunction f(
        [&](const std::vector<double>& x, std::vector<double>& dx) {
            const std::vector<double> resultXX = runSim(compositeXX, x[0], x[1], x[2]);
            const std::vector<double> resultHH = runSim(compositeHH, x[0], x[1], x[2]);
            const std::vector<double> resultXCNOT = runSim(compositeXCNOT, x[0], x[1], x[2]);
            const std::vector<double> resultHCNOT = runSim(compositeHCNOT, x[0], x[1], x[2]);
            double error;
            const auto errorAccum = [&error](const std::vector<double>& in_expected, const std::vector<double>& in_actual){
                for (int i = 0; i < in_expected.size(); ++i)
                {
                    error += std::abs(in_expected[i] - in_actual[i]);
                }
            };
            errorAccum(expectedCompositeXX, resultXX);
            errorAccum(expectedCompositeHH, resultHH);
            errorAccum(expectedCompositeXCNOT, resultXCNOT);
            errorAccum(expectedCompositeHCNOT, resultHCNOT);

            // Returns total/accumulated error
            return error;
        },
        // We have 3 parameters
        3
    );
    
    // Get the XACC nlopt service
    auto optimizer = xacc::getOptimizer("nlopt", { std::make_pair("nlopt-maxeval", 1000), std::make_pair("initial-parameters", std::vector<double>{ 1.25, 0.97, 0.0 }) });
    auto result = optimizer->optimize(f);
   
	std::cout << "Optimal parameters: " << result.second[0] << "; " << result.second[1] << "; " << result.second[2] << "\n";   
    xacc::Finalize();
	return 0;
}