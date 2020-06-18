#include <gtest/gtest.h>
#include "xacc.hpp"
#include "Pulse.hpp"
#include "PulseSystemModel.hpp"
#include "PulseGen.hpp"
#include "CommonGates.hpp"
#include "xacc_service.hpp"

namespace {
    const std::string singleQubitHamiltonianJsonTemplate = R"#(
        {
            "description": "One-qubit transmon qubit Hamiltonian.",
            "h_latex": "",
            "h_str": ["omega0*O0", "g*(SP0 + SM0)||D0", "0.5*delta*O0*(I0-O0)"],
            "osc": {},
            "qub": {
                "0": {{dim}}
            },
            "vars": {
                "omega0": {{omega0}},
                "g": {{g}},
                "delta": {{delta}}
            }
        }
    )#";
    
    std::string createTransmonHamiltonianJson()
    {
        const double omega0 = 5.35*2*M_PI; 
        const double delta = 0.35*2*M_PI; 
        // Driving field coefficient
        const double g = 2.0*M_PI/100;        
        const std::string omega0PlaceHolder = "{{omega0}}";
        const std::string gPlaceHolder = "{{g}}";
        const std::string deltaPlaceHolder = "{{delta}}";
        const std::string dimPlaceHolder = "{{dim}}";

        std::string hamiltonianJson = singleQubitHamiltonianJsonTemplate;
        // Substitute the value to the Hamiltonian
        hamiltonianJson.replace(hamiltonianJson.find(omega0PlaceHolder), omega0PlaceHolder.length(), std::to_string(omega0));
        hamiltonianJson.replace(hamiltonianJson.find(gPlaceHolder), gPlaceHolder.length(), std::to_string(g));
        hamiltonianJson.replace(hamiltonianJson.find(deltaPlaceHolder), deltaPlaceHolder.length(), std::to_string(delta));
        hamiltonianJson.replace(hamiltonianJson.find(dimPlaceHolder), dimPlaceHolder.length(), "2");
        return hamiltonianJson;
    }
}


TEST(DigitalGateTester, testHadamard)
{
    auto systemModel = std::make_shared<QuaC::PulseSystemModel>();
    const bool loadOK = systemModel->loadHamiltonianJson(createTransmonHamiltonianJson());
    assert(loadOK);    
        
    BackendChannelConfigs channelConfigs;
    channelConfigs.dt = 1.0;    
    channelConfigs.loFregs_dChannels.emplace_back(0.0);
    systemModel->setChannelConfigs(channelConfigs);
    
    // We run 10000 shots
    const int NB_SHOTS = 10000;
    auto quaC = xacc::getAccelerator("QuaC", { std::make_pair("system-model", systemModel), std::make_pair("shots", NB_SHOTS)  });    

    auto qubitReg = xacc::qalloc(1);    
    auto xasmCompiler = xacc::getCompiler("xasm");
    auto ir = xasmCompiler->compile(R"(__qpu__ void test1(qbit q) {
      H(q[0]);
      Measure(q[0]);
    })", quaC);

    auto program = ir->getComposite("test1");
    quaC->execute(qubitReg, program);

    const double prob0 = qubitReg->computeMeasurementProbability("0");
    const double prob1 = qubitReg->computeMeasurementProbability("1");
    EXPECT_NEAR(prob0, 0.5, 0.1);
    EXPECT_NEAR(prob1, 0.5, 0.1);
}


TEST(DigitalGateTester, testRx)
{
    auto systemModel = std::make_shared<QuaC::PulseSystemModel>();
    const bool loadOK = systemModel->loadHamiltonianJson(createTransmonHamiltonianJson());
    assert(loadOK);    
        
    BackendChannelConfigs channelConfigs;
    channelConfigs.dt = 1.0;    
    channelConfigs.loFregs_dChannels.emplace_back(0.0);
    systemModel->setChannelConfigs(channelConfigs);
    
    const int NB_SHOTS = 10000;
    auto quaC = xacc::getAccelerator("QuaC", { std::make_pair("system-model", systemModel), std::make_pair("shots", NB_SHOTS)  });    

    auto xasmCompiler = xacc::getCompiler("xasm");
    
    // Rotate Rx(theta) => expected result: cos(theta/2) |0> + sin(theta/2) |1>
    auto ir = xasmCompiler->compile(R"(__qpu__ void test2(qbit q, double t) {
      Rx(q[0], t);
      Measure(q[0]);
    })", quaC);
    auto program = ir->getComposite("test2");
    
    const auto angles = xacc::linspace(-xacc::constants::pi, xacc::constants::pi, 20);
    for (const auto& angle : angles)
    {
        auto qubitReg = xacc::qalloc(1);    
        auto evaled = program->operator()({ angle });
        quaC->execute(qubitReg, evaled);
        const double prob0 = qubitReg->computeMeasurementProbability("0");
        const double prob1 = qubitReg->computeMeasurementProbability("1");
        // cos(theta/2) |0> + sin(theta/2) |1>
        const double expectedProb0 = cos(angle/2)*cos(angle/2);
        const double expectedProb1 = 1.0 - expectedProb0;
        EXPECT_NEAR(prob0, expectedProb0, 0.1);
        EXPECT_NEAR(prob1, expectedProb1, 0.1);
    }
}

TEST(DigitalGateTester, testMultipleGates)
{
    auto systemModel = std::make_shared<QuaC::PulseSystemModel>();
    const bool loadOK = systemModel->loadHamiltonianJson(createTransmonHamiltonianJson());
    assert(loadOK);    
        
    BackendChannelConfigs channelConfigs;
    channelConfigs.dt = 1.0;    
    channelConfigs.loFregs_dChannels.emplace_back(0.0);
    systemModel->setChannelConfigs(channelConfigs);
    
    const int NB_SHOTS = 10000;
    auto quaC = xacc::getAccelerator("QuaC", { std::make_pair("system-model", systemModel), std::make_pair("shots", NB_SHOTS)  });    

    auto xasmCompiler = xacc::getCompiler("xasm");
    
    auto ir = xasmCompiler->compile(R"(__qpu__ void test3(qbit q) {
      H(q[0]);
      Z(q[0]);
      H(q[0]);
      Measure(q[0]);
    })", quaC);
    auto program = ir->getComposite("test3");
    auto qubitReg = xacc::qalloc(1);    
    quaC->execute(qubitReg, program);

    const double prob0 = qubitReg->computeMeasurementProbability("0");
    const double prob1 = qubitReg->computeMeasurementProbability("1");
    // Should be in |1> state
    EXPECT_NEAR(prob0, 0.0, 0.01);
    EXPECT_NEAR(prob1, 1.0, 0.01);
}

TEST(DigitalGateTester, testMultipleQubitGates)
{
    auto systemModel = std::make_shared<QuaC::PulseSystemModel>();
    const bool loadOK = systemModel->loadHamiltonianJson(createTransmonHamiltonianJson());
    assert(loadOK);    
        
    BackendChannelConfigs channelConfigs;
    channelConfigs.dt = 1.0;    
    channelConfigs.loFregs_dChannels.emplace_back(0.0);
    systemModel->setChannelConfigs(channelConfigs);
    
    const int NB_SHOTS = 10000;
    auto quaC = xacc::getAccelerator("QuaC", { std::make_pair("system-model", systemModel), std::make_pair("shots", NB_SHOTS)  });    

    auto xasmCompiler = xacc::getCompiler("xasm");
    
    // This test both the gate decomposition and multi-qubit simulation (gate-based)
    // For example, we have procedures in place to convert a SWAP gate to 3 CNOT gates.
    auto ir = xasmCompiler->compile(R"(__qpu__ void test4(qbit q) {
      X(q[0]);
      Swap(q[0], q[1]);
      Swap(q[1], q[2]);
      Measure(q[0]);
      Measure(q[1]);
      Measure(q[2]);
      Measure(q[3]);
      Measure(q[4]);
    })", quaC);
    auto program = ir->getComposite("test4");
    // We will use 5 qubits and swap the '1' into the middle
    auto qubitReg = xacc::qalloc(5);    
    quaC->execute(qubitReg, program);
    // Qubit in the middle is 1.
    const double probability = qubitReg->computeMeasurementProbability("00100");
    // Probability is 1
    EXPECT_NEAR(probability, 1.0, 0.001);
}

TEST(DigitalGateTester, testMixPulseDigital)
{
    // Similar to the above test: H-Z-H is equiv to an X gate,
    // but use pulses for the Hadamard gate instead.
    
    auto systemModel = std::make_shared<QuaC::PulseSystemModel>();
    const bool loadOK = systemModel->loadHamiltonianJson(createTransmonHamiltonianJson());
    assert(loadOK);   
    const int total_samples = 100; 
    BackendChannelConfigs channelConfigs;
    channelConfigs.dt = 1.0;    
    channelConfigs.loFregs_dChannels.emplace_back(5.35);
    
    // For this system: driving strength is set to
    // const double g = 2.0*M_PI/100;  
    // i.e. 100 sample = 2pi rotation
    // we want pi/2 => square pulse of 25 samples
    channelConfigs.addOrReplacePulse("square", QuaC::SquarePulse(25));     
    systemModel->setChannelConfigs(channelConfigs);
    
    // We run 10000 shots
    const int NB_SHOTS = 10000;
    auto quaC = xacc::getAccelerator("QuaC", { std::make_pair("system-model", systemModel), std::make_pair("shots", NB_SHOTS)  });    

    auto qubitReg = xacc::qalloc(1);    
    auto provider = xacc::getIRProvider("quantum");
    auto compositeInst = provider->createComposite("test_hadamard_square_pulse");
    
    // First Hadamard (Pulse)
    {
        // We need to frame change in order to create a Hadarmard
        {
            auto fcInst = std::make_shared<xacc::quantum::Pulse>("fc", "d0");
            xacc::InstructionParameter fcParameter(-M_PI/2);
            fcInst->setParameter(0, fcParameter);
            fcInst->setBits({0});
            compositeInst->addInstruction(fcInst);
        }
        {
            // Add X(pi/2) pulse (Gaussian with optimal pulse width)
            auto pulseInst = std::make_shared<xacc::quantum::Pulse>("square", "d0");
            pulseInst->setBits({0});
            // Add the Gaussian pulse
            compositeInst->addInstruction(pulseInst);
        }
        {
            // Another frame-change
            auto fcInst = std::make_shared<xacc::quantum::Pulse>("fc", "d0");
            xacc::InstructionParameter fcParameter(-M_PI/2);
            fcInst->setParameter(0, fcParameter);
            fcInst->setBits({0});
            compositeInst->addInstruction(fcInst);
        }
    }
        
    // We should have applied a Hadamard gate by *pulse*
    // Now apply a digital Z gate:
    {
        auto zGate = std::make_shared<xacc::quantum::Z>(0);
        compositeInst->addInstruction(zGate);
    }
    
    // Next: another Hadamard gate: by pulse
    {
        {
            auto fcInst = std::make_shared<xacc::quantum::Pulse>("fc", "d0");
            xacc::InstructionParameter fcParameter(-M_PI/2);
            fcInst->setParameter(0, fcParameter);
            fcInst->setBits({0});
            compositeInst->addInstruction(fcInst);
        }
        {
            // Add X(pi/2) pulse (Gaussian with optimal pulse width)
            auto pulseInst = std::make_shared<xacc::quantum::Pulse>("square", "d0");
            pulseInst->setBits({0});
            // Add the Gaussian pulse
            compositeInst->addInstruction(pulseInst);
        }
        {
            // Another frame-change
            auto fcInst = std::make_shared<xacc::quantum::Pulse>("fc", "d0");
            xacc::InstructionParameter fcParameter(-M_PI/2);
            fcInst->setParameter(0, fcParameter);
            fcInst->setBits({0});
            compositeInst->addInstruction(fcInst);
        }
    }
   
    // Measure qubit after the pulse
    auto meas = std::make_shared<xacc::quantum::Measure>(0);
    compositeInst->addInstruction(meas);

    // Run the Pulse simulation with the Hamiltonian provided
    quaC->execute(qubitReg, compositeInst);
    
    const double prob0 = qubitReg->computeMeasurementProbability("0");
    const double prob1 = qubitReg->computeMeasurementProbability("1");
    // Should be in |1> state
    EXPECT_NEAR(prob0, 0.0, 0.1);
    EXPECT_NEAR(prob1, 1.0, 0.1);  
}

TEST(DigitalGateTester, testQuantumProcessTomography)
{
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

    auto qpt = xacc::getService<xacc::Algorithm>("qpt");
    // First: use Qpp backend to get the *reference* data
    auto qppReg = xacc::qalloc(1);    
    // Use QPP accelerator   
    auto acc = xacc::getAccelerator("qpp", {std::make_pair("shots", 2048)});
    auto compiler = xacc::getCompiler("xasm");
    auto ir = compiler->compile(R"(__qpu__ void f(qbit q) {
        X(q[0]);
    })", nullptr);
    auto qppCompositeInstr = ir->getComposite("f");
    qpt->initialize({ std::make_pair("circuit", qppCompositeInstr), std::make_pair("accelerator", acc)});
    qpt->execute(qppReg);
    qppReg->print();
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
    auto quaC = xacc::getAccelerator("QuaC", { std::make_pair("system-model", systemModel), std::make_pair("shots", NB_SHOTS)  });    
    // Perform a QPT with QuaC
    qpt->initialize({ std::make_pair("circuit", compositeInst), std::make_pair("accelerator", quaC), std::make_pair("optimize-circuit", false)});
    qpt->execute(qubitReg);    
    qubitReg->print();

    // QuaC result data:
    const auto quac_chi_real_vec = (*qubitReg)["chi-real"].as<std::vector<double>>();
    const auto quac_chi_imag_vec = (*qubitReg)["chi-imag"].as<std::vector<double>>(); 
    
    EXPECT_EQ(qpp_chi_real_vec.size(), qpp_chi_imag_vec.size());
    EXPECT_EQ(qpp_chi_real_vec.size(), quac_chi_real_vec.size());
    EXPECT_EQ(qpp_chi_imag_vec.size(), quac_chi_imag_vec.size());

    for (int i = 0; i < qpp_chi_real_vec.size(); ++i)
    {
        // Check that all elements in the chi vector are matched.
        EXPECT_NEAR(qpp_chi_real_vec[i], quac_chi_real_vec[i], 0.1);
        EXPECT_NEAR(qpp_chi_imag_vec[i], quac_chi_imag_vec[i], 0.1);
    }
}

// Test digital simulation (Deuteuron)
TEST(DigitalGateTester, testDeuteuron)
{
    const int NB_SHOTS = 10000;
    auto systemModel = std::make_shared<QuaC::PulseSystemModel>();
    const bool loadOK = systemModel->loadHamiltonianJson(createTransmonHamiltonianJson());
    assert(loadOK);    
    BackendChannelConfigs channelConfigs;
    channelConfigs.dt = 1.0;    
    channelConfigs.loFregs_dChannels.emplace_back(0.0);
    systemModel->setChannelConfigs(channelConfigs);
    // Use QuaC
    auto accelerator = xacc::getAccelerator("QuaC", { std::make_pair("system-model", systemModel), std::make_pair("shots", NB_SHOTS)  });    
    auto xasmCompiler = xacc::getCompiler("xasm");
    auto ir = xasmCompiler->compile(R"(__qpu__ void ansatz(qbit q, double t) {
      X(q[0]);
      Ry(q[1], t);
      CX(q[1], q[0]);
      H(q[0]);
      H(q[1]);
      Measure(q[0]);
      Measure(q[1]);
    })", accelerator);

    auto program = ir->getComposite("ansatz");
    // Expected results from deuteron_2qbit_xasm_X0X1
    const std::vector<double> expectedResults {
        0.0,
        -0.324699,
        -0.614213,
        -0.837166,
        -0.9694,
        -0.996584,
        -0.915773,
        -0.735724,
        -0.475947,
        -0.164595,
        0.164595,
        0.475947,
        0.735724,
        0.915773,
        0.996584,
        0.9694,
        0.837166,
        0.614213,
        0.324699,
        0.0
    };

    const auto angles = xacc::linspace(-xacc::constants::pi, xacc::constants::pi, 20);
    for (size_t i = 0; i < angles.size(); ++i)
    {
        auto buffer = xacc::qalloc(2);
        auto evaled = program->operator()({ angles[i] });
        accelerator->execute(buffer, evaled);
        // Error bar limit: 0.1
        EXPECT_NEAR(buffer->getExpectationValueZ(), expectedResults[i], 0.1);
    }
}

int main(int argc, char **argv) 
{
  xacc::Initialize(argc, argv);
  ::testing::InitGoogleTest(&argc, argv);
  auto ret = RUN_ALL_TESTS();
  xacc::Finalize();
  return ret;
}