#include "xacc.hpp"
#include <cassert>
#include "Pulse.hpp"
#include "xacc_service.hpp"
#include "PulseSystemModel.hpp"
#include "PulseGen.hpp"
#include <chrono>

// This demonstrates the procedure to calibrate a X(pi/2) pulse
// which is the key pulse to construct single-qubit gates.
int main (int argc, char** argv) {
	// Initialize the XACC Framework
	xacc::Initialize(argc, argv);
    // Expected chi matrix
    std::vector<double> chi_exp_real_vec;
    std::vector<double> chi_exp_imag_vec;

    {
        auto systemModel = std::make_shared<QuaC::PulseSystemModel>();
    
        const std::string hamiltonianJson = R"#(
            {
                "description": "Two-qubit Hamiltonian",
                "h_str": ["_SUM[i,0,1,wq{i}*O{i}]", "_SUM[i,0,1,delta{i}*O{i}*(O{i}-I{i})]", "_SUM[i,0,1,omegad{i}*X{i}||D{i}]", "omegad1*X0||U0", "omegad0*X1||U1", "jq0q1*Sp0*Sm1", "jq0q1*Sm0*Sp1"],
                "osc": {},
                "qub": {
                    "0": 2,
                    "1": 2
                },
                "vars": {
                    "wq0": 30.518812656662774, 
                    "wq1": 31.238229295532093,
                    "delta0": -2.011875935,
                    "delta1": -2.008734343,
                    "omegad0": -1.703999855,
                    "omegad1": -1.703999855,
                    "jq0q1": 0.011749557 
                }
            }
        )#";

        // Load Hamiltonian
        systemModel->loadHamiltonianJson(hamiltonianJson);

        std::vector<double> d_loFreqs { 4.857, 4.972 };
        std::vector<double> u_loFreqs { 4.972, 4.857 };

        BackendChannelConfigs channelConfigs;
        channelConfigs.dt = 0.222;
        channelConfigs.loFregs_dChannels = d_loFreqs;
        channelConfigs.loFregs_uChannels = u_loFreqs;
        systemModel->setChannelConfigs(channelConfigs);    

        // Compute the expected chi matrix (Rx(pi/2))
        auto acc = xacc::getAccelerator("QuaC", { std::make_pair("system-model", systemModel), std::make_pair("shots", 10000) });    
        auto buffer = xacc::qalloc(2);
        auto compiler = xacc::getCompiler("xasm");
        auto ir = compiler->compile(R"(__qpu__ void rx(qbit q) {
            Rx(q[0], pi/2);
        })", nullptr);
        auto rx = ir->getComposite("rx");
        auto qpt = xacc::getService<xacc::Algorithm>("qpt");
        qpt->initialize({ 
            std::make_pair("circuit", rx), 
            std::make_pair("accelerator", acc), 
            std::make_pair("optimize-circuit", false)
        });
        qpt->execute(buffer);

        chi_exp_real_vec = (*buffer)["chi-real"].as<std::vector<double>>();
        chi_exp_imag_vec = (*buffer)["chi-imag"].as<std::vector<double>>();

        std::cout << "Expeted Chi Matrix:\n";
        for (int i = 0; i < chi_exp_real_vec.size(); ++i)
        {
            std::cout << "(" << chi_exp_real_vec[i] << ", " << chi_exp_imag_vec [i] << "), ";
        }
        std::cout << "\n\n";
    }

    // Calibrate Rx(pi/2)
    xacc::OptFunction qptOptFunc([&](const std::vector<double>& x, std::vector<double>& dx) {
        const std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
        
        auto systemModel = std::make_shared<QuaC::PulseSystemModel>();
    
        const std::string hamiltonianJson = R"#(
            {
                "description": "Two-qubit Hamiltonian",
                "h_str": ["_SUM[i,0,1,wq{i}*O{i}]", "_SUM[i,0,1,delta{i}*O{i}*(O{i}-I{i})]", "_SUM[i,0,1,omegad{i}*X{i}||D{i}]", "omegad1*X0||U0", "omegad0*X1||U1", "jq0q1*Sp0*Sm1", "jq0q1*Sm0*Sp1"],
                "osc": {},
                "qub": {
                    "0": 2,
                    "1": 2
                },
                "vars": {
                    "wq0": 30.518812656662774, 
                    "wq1": 31.238229295532093,
                    "delta0": -2.011875935,
                    "delta1": -2.008734343,
                    "omegad0": -1.703999855,
                    "omegad1": -1.703999855,
                    "jq0q1": 0.011749557 
                }
            }
        )#";

        // Load Hamiltonian
        systemModel->loadHamiltonianJson(hamiltonianJson);

        std::vector<double> d_loFreqs { 4.857, 4.972 };
        std::vector<double> u_loFreqs { 4.972, 4.857 };

        const int nSamples = 121;
        BackendChannelConfigs channelConfigs;
        channelConfigs.dt = 0.222;
        channelConfigs.loFregs_dChannels = d_loFreqs;
        channelConfigs.loFregs_uChannels = u_loFreqs;

        
        const double ampl = x[0];
        if (ampl < 0.05 || ampl > 0.2)
        {
            return 1.0;
        }

        const double beta = 1.0;
        const std::string pulseName = "drag_" + std::to_string(ampl) + "_" + std::to_string(beta);
        channelConfigs.addOrReplacePulse(pulseName, QuaC::Drag(nSamples, ampl, 20, beta));
        systemModel->setChannelConfigs(channelConfigs);    

        auto provider = xacc::getIRProvider("quantum");
        auto composite = provider->createComposite("test_pulse");
        const std::string channel = "d0";
        auto pulse = std::make_shared<xacc::quantum::Pulse>(pulseName, channel);
        pulse->setBits({0});

        // Frame-change (by an angle to be optimized for)   
        auto fcInst = std::make_shared<xacc::quantum::Pulse>("fc", "d0");
        xacc::InstructionParameter fcParameter(x[1]);
        fcInst->setParameter(0, fcParameter);
        
        composite->addInstruction(fcInst);
        composite->addInstruction(pulse);
        auto qubitReg = xacc::qalloc(2);
        const int NB_SHOTS = 10000;
        auto qpu = xacc::getAccelerator("QuaC", { std::make_pair("system-model", systemModel), std::make_pair("shots", NB_SHOTS)  });    

        // Get Quantum Process Tomography Algo
        auto qpt = xacc::getAlgorithm("qpt");
        qpt->initialize({
            std::make_pair("circuit", composite), 
            std::make_pair("accelerator", qpu), 
            std::make_pair("optimize-circuit", false)
        });
        qpt->execute(qubitReg);    


        auto chi_real_vec = (*qubitReg)["chi-real"].as<std::vector<double>>();
        auto chi_imag_vec = (*qubitReg)["chi-imag"].as<std::vector<double>>();

        std::cout << "Result Chi Matrix:\n";
        for (int i = 0; i < chi_real_vec.size(); ++i)
        {
            std::cout << "(" << chi_real_vec[i] << ", " << chi_imag_vec [i] << "), ";
        }

        double normDiff = 0.0;
        for (int i = 0; i < chi_real_vec.size(); ++i)
        {
            const std::complex<double> actual { chi_real_vec[i],  chi_imag_vec[i] };
            const std::complex<double> theory { chi_exp_real_vec[i],  chi_exp_imag_vec[i] };
            normDiff += std::norm(actual - theory);
        }


        // Calculates the fidelity based on QPT result
        xacc::HeterogeneousMap chiReference { 
            std::make_pair("chi-theoretical-real", chi_exp_real_vec), 
            std::make_pair("chi-theoretical-imag", chi_exp_imag_vec) };

        const double fidelityResult = qpt->calculate("fidelity", qubitReg, chiReference);

        // Print out for tracking
        const std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();    
        std::cout << "=====================================================\n";
        std::cout << "Elapsed time = " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << " [secs]\n";
        std::cout << "\n";
        std::cout << "Fidelity result: " << fidelityResult << "\n";
        std::cout << "Norm diff: " << normDiff << "\n";

        std::cout << "=====================================================\n";
        
        return normDiff;
    }, 
    2);
    
    // Using nl-opt optimizer to optimize 
    auto optimizer = xacc::getOptimizer("nlopt", { std::make_pair("nlopt-maxeval", 50) });

    auto result = optimizer->optimize(qptOptFunc);
    std::cout << "Final fidelity error: " << result.first << "\n";
    std::cout << "Optimal params: ";
    for (const auto& par : result.second)
    {
        std::cout << par << ", ";
    }

    xacc::Finalize();
	return 0;
}