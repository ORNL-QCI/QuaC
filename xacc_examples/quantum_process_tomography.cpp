#include "xacc.hpp"
#include "xacc_service.hpp"
#include "Pulse.hpp"
#include "PulseSystemModel.hpp"
#include "PulseGen.hpp"

namespace {
    const std::string hamiltonianJson =  R"(
        {
            "description": "One-qubit Hamiltonian.",
            "h_latex": "",
            "h_str": ["0.5*epsilon*Z0", "0.5*delta*X0||D0"],
            "osc": {},
            "qub": {
                "0": 2
            },
            "vars": {
                "epsilon": 0.0,
                "delta": 6.2831853
            }
        }
    )";
}

// In this example, we'll use the XACC Quantum Process Tomography plugin
// to compare a quantum gate process created by pulse vs. that of an ideal gate.
// In this example, we demonstrate QPT for a simple X gate (Pi pulse).
int main (int argc, char** argv) {
	// Initialize the XACC Framework
	xacc::Initialize(argc, argv);
    auto qpt = xacc::getService<xacc::Algorithm>("qpt");

    // First: use Qpp backend to get the *reference* data
    auto qppReg = xacc::qalloc(1);    
    // Use QPP accelerator   
    auto acc = xacc::getAccelerator("qpp", {std::make_pair("shots", 1024)});
    auto compiler = xacc::getCompiler("xasm");
    auto ir = compiler->compile(R"(__qpu__ void f(qbit q) {
        X(q[0]);
    })", nullptr);
    auto qppCompositeInstr = ir->getComposite("f");
    qpt->initialize({ std::make_pair("circuit", qppCompositeInstr), std::make_pair("accelerator", acc)});
    qpt->execute(qppReg);
    // Reference data:
    const auto qpp_chi_real_vec = (*qppReg)["chi-real"].as<std::vector<double>>();
    const auto qpp_chi_imag_vec = (*qppReg)["chi-imag"].as<std::vector<double>>(); 

    // Now, let's use QuaC to perform a Pi pulse (X gate)
    // Load Hamiltonian model:
    auto systemModel = std::make_shared<QuaC::PulseSystemModel>();
    systemModel->loadHamiltonianJson(hamiltonianJson);
    BackendChannelConfigs channelConfigs;
    channelConfigs.dt = 0.005;
    channelConfigs.loFregs_dChannels.emplace_back(0.0);
    const size_t nbSamples = 100;

    // Just use a simple Pi pulse
    channelConfigs.addOrReplacePulse("square", QuaC::SquarePulse(nbSamples));
    systemModel->setChannelConfigs(channelConfigs);    
    auto qubitReg = xacc::qalloc(1);    
    auto provider = xacc::getIRProvider("quantum");
    auto compositeInst = provider->createComposite("pulse_qpt");
    // X Pi pulse
    auto pulseInst = std::make_shared<xacc::quantum::Pulse>("square", "d0");
    pulseInst->setBits({0});
    compositeInst->addInstruction(pulseInst);

    // We run 10000 shots
    const int NB_SHOTS = 10000;
    // Note: since this circuit contains Pulse instructions (constructed via IR),
    // we don't want to perform circuit optimization on the QPT test circuits, hence set optimize-circuit flag to false.
    // (e.g. circuit optimization is based on graph which may not be able to handle circuits that have both pulse and gate instructions)
    auto quaC = xacc::getAccelerator("QuaC", { std::make_pair("system-model", systemModel), std::make_pair("shots", NB_SHOTS) });    
    // Perform a QPT with QuaC
    qpt->initialize({ std::make_pair("circuit", compositeInst), std::make_pair("accelerator", quaC), std::make_pair("optimize-circuit", false)});
    qpt->execute(qubitReg);    

    // Now, let's calculate the Fidelity b/w QuaC and reference (Qpp)
    xacc::HeterogeneousMap chiReference { std::make_pair("chi-theoretical-real", qpp_chi_real_vec), std::make_pair("chi-theoretical-imag", qpp_chi_imag_vec) };
    const double result = qpt->calculate("fidelity", qubitReg, chiReference);

    std::cout << "Fidelity result: " << result << "\n";
    // Finalize the XACC Framework
	xacc::Finalize();

	return 0;
}