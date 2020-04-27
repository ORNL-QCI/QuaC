#include "xacc.hpp"
#include <cassert>
#include "Pulse.hpp"
#include "xacc_service.hpp"
#include "PulseSystemModel.hpp"
#include "PulseGen.hpp"
#include <chrono>

namespace {
struct u3_param
{
    double theta;
    double phi;
    double lambda;
    u3_param (double in_theta, double in_phi, double in_lambda):
        theta(in_theta), 
        phi(in_phi), 
        lambda(in_lambda)
    {}

    std::string toString() const 
    {
        return "U3(" + std::to_string(theta) + ", " + std::to_string(phi) + ", " + std::to_string(lambda) + ")";
    }
};
}


int main (int argc, char** argv) {
	// Initialize the XACC Framework
	xacc::Initialize(argc, argv);

    xacc::OptFunction qptOptFunc([&](const std::vector<double>& x, std::vector<double>& dx) {
        // Time the execution 
        const std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
        
        std::vector<u3_param> params;
        for (int i = 0; i < x.size(); i += 3)
        {
            params.emplace_back(x[i], x[i+1], x[i+2]);
        }
        // Expects 4 U3 gates
        assert(params.size() == 4);
        
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
        //  CR drive length = 848*dt
        const int nSamples = 848;
        BackendChannelConfigs channelConfigs;
        channelConfigs.dt = 0.222;
        channelConfigs.loFregs_dChannels = d_loFreqs;
        channelConfigs.loFregs_uChannels = u_loFreqs;

        // Amplitude: see IPython notebook to see how we calibrate this amplitude.
        // CR_01: 0.47307243770437246 
        const double A = 0.47307243770437246;
        const std::string pulseName = "square" + std::to_string(A);
        channelConfigs.addOrReplacePulse(pulseName, QuaC::GaussianSquare(nSamples, A, 32, 720));
        systemModel->setChannelConfigs(channelConfigs);    

        auto provider = xacc::getIRProvider("quantum");
        auto composite = provider->createComposite("test_pulse");
        const std::string channel = "u0";
        auto pulse = std::make_shared<xacc::quantum::Pulse>(pulseName, channel);
        pulse->setBits({0, 1});

        // Pre- and Post- U3 gates
        auto q0_pre = provider->createInstruction("U", {0}, { params[0].theta, params[0].phi, params[0].lambda });
        auto q1_pre = provider->createInstruction("U", {1}, { params[1].theta, params[1].phi, params[1].lambda });
        auto q0_post = provider->createInstruction("U", {0}, { params[2].theta, params[2].phi, params[2].lambda });
        auto q1_post = provider->createInstruction("U", {1}, { params[3].theta, params[3].phi, params[3].lambda });
        
        // Construct the circuit for QPT
        // Pre-pulse U3 gates
        composite->addInstruction(q0_pre);
        composite->addInstruction(q1_pre);
        // CR pulse
        composite->addInstruction(pulse);
        // Post-pulse U3 gates
        composite->addInstruction(q0_post);
        composite->addInstruction(q1_post);

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
        // Target Chi matrix of CNOT gate
        const std::vector<double> true_cx_chi_real {
        1., 0., 0., 1., 1., 0., 0., -1., 0., 0., 0., 0., 0., 0., 0.,
        0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
        0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0.,
        1., 1., 0., 0., -1., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 1., 1.,
        0., 0., -1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
        0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
        0., 0., 0., 0., 0., 0., 0., -1., 0., 0., -1., -1., 0., 0., 1., 0., 0.,
        0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
        0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
        0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
        0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
        0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
        0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
        0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
        0., 0., 0., 0., 0., 0., 0., 0.
        };
        assert(true_cx_chi_real.size() == 256);
        // Chi imag: all zero's
        const std::vector<double> true_cx_chi_imag(256, 0.0);

        // Calculates the fidelity based on QPT result
        xacc::HeterogeneousMap chiReference { std::make_pair("chi-theoretical-real", true_cx_chi_real), std::make_pair("chi-theoretical-imag", true_cx_chi_imag) };
        const double fidelityResult = qpt->calculate("fidelity", qubitReg, chiReference);

        // Print out for tracking
        const std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();    
        std::cout << "=====================================================\n";
        std::cout << "Elapsed time = " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << " [secs]\n";
        for (const auto& u3Par : params)
        {
            std::cout << u3Par.toString() << ", ";
        }
        std::cout << "\n";
        std::cout << "Fidelity result: " << fidelityResult << "\n";
        std::cout << "=====================================================\n";

        // Target fidelity is 1.0 => cost function = abs(1.0 - fidelity)
        return std::fabs(1.0 - fidelityResult);
    }, 
    // We have 12 params
    12);
    
    // Using nl-opt optimizer to optimize these 12 params for U3 gates:
    auto optimizer = xacc::getOptimizer("nlopt", { std::make_pair("nlopt-maxeval", 1000) });

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