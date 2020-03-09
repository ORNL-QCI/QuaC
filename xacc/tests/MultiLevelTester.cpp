#include <gtest/gtest.h>
#include "xacc.hpp"
#include "Pulse.hpp"
#include "PulseSystemModel.hpp"
#include "PulseGen.hpp"
#include "CommonGates.hpp"

namespace {
    // Three-level single-qubit Hamiltonian
    // This is the same as a regular 2-level Hamiltonian but Pauli X and Z operators
    // are converted to ladder operators.
    const std::string singleQubitNonPauliTmpl = R"(
        {
            "description": "One-qubit Hamiltonian.",
            "h_latex": "",
            "h_str": ["-0.5*omega0*I0", "omega0*O0", "omegaa*(SM0+SP0)||D0"],
            "osc": {},
            "qub": {
                "0": 3
            },
            "vars": {
                "omega0": {{omega0}},
                "omegaa": {{omegaa}}
            }
        }
    )";

    // Qutrit (dim=3) Pauli form:
    // QuaC backend can now handle conversion from Pauli ops to ladder ops for
    // higher dimensional qubit systems.  
    const std::string singleQubitPauliTmpl = R"(
        {
            "description": "One-qubit Hamiltonian.",
            "h_latex": "",
            "h_str": ["-0.5*omega0*Z0", "omegaa*X0||D0"],
            "osc": {},
            "qub": {
                "0": 3
            },
            "vars": {
                "omega0": {{omega0}},
                "omegaa": {{omegaa}}
            }
        }
    )";

    std::string createHamiltonianJson1Q(double in_omega0, double in_omegaA, bool in_PauliForm)
    {
        const std::string omega0PlaceHolder = "{{omega0}}";
        const std::string omegaaPlaceHolder = "{{omegaa}}";
        std::string hamiltonianJson = in_PauliForm ? singleQubitPauliTmpl : singleQubitNonPauliTmpl;
        hamiltonianJson.replace(hamiltonianJson.find(omega0PlaceHolder), omega0PlaceHolder.length(), std::to_string(in_omega0));
        hamiltonianJson.replace(hamiltonianJson.find(omegaaPlaceHolder), omegaaPlaceHolder.length(), std::to_string(in_omegaA));
        return hamiltonianJson;
    }
}


TEST(MultiLevelTester, testThreeLevel)
{
    const auto runTest = [](bool pauliForm){
        auto systemModel = std::make_shared<QuaC::PulseSystemModel>();
        // Params:
        const int total_samples = 100;
        const double omega_0 = 2 * M_PI;
        const double omega_d0 = omega_0;
        
        // Max |1> population: 
        // Expected: Pop(0) = 4/9; Pop(1) = 1/3; Pop(2) = 2/9 
        const double omega_a = M_PI / (std::sqrt(3.0) * total_samples);

        const std::string hamiltonianJson = createHamiltonianJson1Q(omega_0, omega_a, pauliForm);

        const bool loadOK = systemModel->loadHamiltonianJson(hamiltonianJson);
        EXPECT_TRUE(loadOK);
        BackendChannelConfigs channelConfigs;
        
        channelConfigs.dt = 1.0;
        // LO freq: this is non-zero, hence the drive signal will be modulated. 
        channelConfigs.loFregs_dChannels.emplace_back(omega_d0/(2*M_PI));
        
        // Add the square pulse
        channelConfigs.addOrReplacePulse("square", QuaC::SquarePulse(total_samples));    
        systemModel->setChannelConfigs(channelConfigs);
        const int NB_SHOTS = 1000;
        auto quaC = xacc::getAccelerator("QuaC", { std::make_pair("system-model", systemModel), std::make_pair("shots", NB_SHOTS) });    

        auto qubitReg = xacc::qalloc(1);    

        auto provider = xacc::getIRProvider("quantum");
        auto compositeInst = provider->createComposite("test_pulse");
        auto pulseInst = std::make_shared<xacc::quantum::Pulse>("square", "d0");
        compositeInst->addInstruction(pulseInst);
        
        // Run the Pulse simulation with the Hamiltonian provided
        quaC->execute(qubitReg, compositeInst);
        qubitReg->print();
        const auto finalPopulations = qubitReg->getInformation("DensityMatrixDiags").as<std::vector<double>>();
        // We expect to have 3 numbers
        EXPECT_EQ(finalPopulations.size(), 3);
        // Expected: Pop(0) = 4/9; Pop(1) = 1/3; Pop(2) = 2/9 
        EXPECT_NEAR(finalPopulations[0], 4.0/9.0, 0.01);
        EXPECT_NEAR(finalPopulations[1], 1.0/3.0, 0.01);
        EXPECT_NEAR(finalPopulations[2], 2.0/9.0, 0.01);
    };

    runTest(false);
    runTest(true);
}

// Transmon qubit Hamiltonian (three-level)
// in the ladder operator form.
TEST(MultiLevelTester, transmonQubitHamiltonian)
{
    
    // Three-level single-qubit Hamiltonian
    // Ref: Phys. Rev. A 96, 022330 
    // Equation (6), params in Sec III, first paragraph
    const std::string transmonHamiltonian = R"#(
        {
            "description": "One-qubit Hamiltonian.",
            "h_latex": "",
            "h_str": ["(w - 0.5*alpha)*O0", "0.5*alpha*O0*O0", "O*(SM0 + SP0)||D0"],
            "osc": {},
            "qub": {
                "0": 3
            },
            "vars": {
                "w": 31.63772297724,
                "alpha": -1.47969,
                "O": 0.0314
            }
        }
    )#";
    
    auto systemModel = std::make_shared<QuaC::PulseSystemModel>();
    const bool loadOk = systemModel->loadHamiltonianJson(transmonHamiltonian);
    EXPECT_TRUE(loadOk);
    
    BackendChannelConfigs channelConfigs;
    channelConfigs.dt = 1.0;
    channelConfigs.loFregs_dChannels.emplace_back(5.0353);
    // 100 => PI pulse
    const size_t nbSamples = 100;
    // Constant drive
    channelConfigs.addOrReplacePulse("square", QuaC::SquarePulse(nbSamples));
    
    systemModel->setChannelConfigs(channelConfigs);
    auto quaC = xacc::getAccelerator("QuaC", { std::make_pair("system-model", systemModel) });    
    auto qubitReg = xacc::qalloc(1);    

    auto provider = xacc::getIRProvider("quantum");
    auto compositeInst = provider->createComposite("test_pulse");
    auto pulseInstD0 = std::make_shared<xacc::quantum::Pulse>("square", "d0");

    compositeInst->addInstruction(pulseInstD0);

    // Run the Pulse simulation with the Hamiltonian provided
    quaC->execute(qubitReg, compositeInst);
    qubitReg->print();
    const auto finalPopulations = qubitReg->getInformation("DensityMatrixDiags").as<std::vector<double>>();
    // We expect to have 3 numbers
    EXPECT_EQ(finalPopulations.size(), 3);
    // PI pulse
    EXPECT_NEAR(finalPopulations[0], 0.0, 0.01);
    EXPECT_NEAR(finalPopulations[1], 1.0, 0.01);
    // Not so much leakage is expected
    EXPECT_NEAR(finalPopulations[2], 0.0, 0.01);
}

// Test cross-resonance b/w two transmon qubits
// modeled as three-level system 
// Note: use target Resonator's Reference Frame + RWA
// (Eugene's Qutip Ipython notebook, second section)
TEST(MultiLevelTester, crossResonance)
{
    const auto runTest = [](std::shared_ptr<xacc::AcceleratorBuffer>& io_buffer, bool in_controlState, int in_nbSamples){
        // Note: because we are in the target qubit (Q1) frame,
        // the self-energy term of Q0 is (w_0 - w_1)
        const std::string rotatingFrameHamiltonian =  R"#(
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
        
        auto systemModel = std::make_shared<QuaC::PulseSystemModel>();
        const bool loadOk = systemModel->loadHamiltonianJson(rotatingFrameHamiltonian);
        EXPECT_TRUE(loadOk);
        
        BackendChannelConfigs channelConfigs;
        channelConfigs.dt = 1.0;
        // LO freq = 0.0 (rotating frame)
        channelConfigs.loFregs_dChannels.emplace_back(0.0);
        const std::string pulseName = "square" + std::to_string(in_nbSamples);
        // Constant drive
        channelConfigs.addOrReplacePulse(pulseName, QuaC::SquarePulse(in_nbSamples));
        
        systemModel->setChannelConfigs(channelConfigs);
        // If the control qubit is set, assign the initial state to 1
        if (in_controlState)
        {
            systemModel->setQubitInitialPopulation(0, 1.0);
        }

        auto quaC = xacc::getAccelerator("QuaC", { std::make_pair("system-model", systemModel) });    
        auto provider = xacc::getIRProvider("quantum");
        auto compositeInst = provider->createComposite("test_pulse" + pulseName);
        auto pulseInstD0 = std::make_shared<xacc::quantum::Pulse>(pulseName, "d0");
        compositeInst->addInstruction(pulseInstD0);
        // Run the Pulse simulation with the Hamiltonian provided
        quaC->execute(io_buffer, compositeInst);
    };
    
    // ================== Control = 0 ======================================================================
    // Expectation Z of target qubit oscillate at slow speed (period ~ 4000)
    // Test 1: control qubit is 0: we drive to t = 2000, target qubit expect-Z -> -1 (see Notebook)
    {
        auto qubitReg = xacc::qalloc(2);    
        runTest(qubitReg, false, 2000);
        qubitReg->print();
        const auto densityMatrixDiags = qubitReg->getInformation("DensityMatrixDiags").as<std::vector<double>>();
        // The density matrix size is 9 x 9 
        EXPECT_EQ(densityMatrixDiags.size(), 9);
        // We have 2 qubit systems
        const auto finalPopulations = qubitReg->getInformation("<O>").as<std::vector<double>>();
        EXPECT_EQ(finalPopulations.size(), 2);
        // Qubit 1 is in ~|1> state @ t = 2000
        EXPECT_NEAR(finalPopulations[1], 1.0, 0.1);
        // Qubit 0 (control) is near 0.0 
        // (note: there are some oscillation in the control qubit as can be seen in the notebook)    
        EXPECT_NEAR(finalPopulations[0], 0.0, 0.25);
    }
   
    // Test 2: control qubit is 0: we drive to t = 3900, target qubit expect-Z -> 1 (see Notebook)
    // This completes the test for the case where control qubit is 0
    {
        auto qubitReg = xacc::qalloc(2);    
        runTest(qubitReg, false, 3900);
        qubitReg->print();
        const auto densityMatrixDiags = qubitReg->getInformation("DensityMatrixDiags").as<std::vector<double>>();
        // The density matrix size is 9 x 9 
        EXPECT_EQ(densityMatrixDiags.size(), 9);
        // We have 2 qubit systems
        const auto finalPopulations = qubitReg->getInformation("<O>").as<std::vector<double>>();
        EXPECT_EQ(finalPopulations.size(), 2);
        // Qubit 1 is in ~|0> state @ t = 3900
        EXPECT_NEAR(finalPopulations[1], 0.0, 0.1);
        // Qubit 0 (control) is near 0.0 
        // (note: there are some oscillation in the control qubit as can be seen in the notebook)    
        EXPECT_NEAR(finalPopulations[0], 0.0, 0.25);
    }
    // ===================================================================================================
    // ============================ Control = 1 ===========================================================
    // Expectation Z of target qubit oscillate at high speed (period ~ 1700)
    // Test 1: control qubit is 1: we drive to t = 900, target qubit expect-Z -> -1 (see Notebook)
    {
        auto qubitReg = xacc::qalloc(2);    
        runTest(qubitReg, true, 900);
        qubitReg->print();
        const auto densityMatrixDiags = qubitReg->getInformation("DensityMatrixDiags").as<std::vector<double>>();
        // The density matrix size is 9 x 9 
        EXPECT_EQ(densityMatrixDiags.size(), 9);
        // We have 2 qubit systems
        const auto finalPopulations = qubitReg->getInformation("<O>").as<std::vector<double>>();
        EXPECT_EQ(finalPopulations.size(), 2);
        // Qubit 1 is in ~|1> state @ t = 900
        EXPECT_NEAR(finalPopulations[1], 1.0, 0.1);
        // Qubit 0 (control) is near 1.0 
        // (note: there are some oscillation in the control qubit as can be seen in the notebook)    
        EXPECT_NEAR(finalPopulations[0], 1.0, 0.25);
    }
   
    // Test 2: control qubit is 1: we drive to t = 1700, target qubit expect-Z -> 1 (see Notebook)
    // This completes the test for the case where control qubit is 1
    {
        auto qubitReg = xacc::qalloc(2);    
        runTest(qubitReg, true, 1700);
        qubitReg->print();
        const auto densityMatrixDiags = qubitReg->getInformation("DensityMatrixDiags").as<std::vector<double>>();
        // The density matrix size is 9 x 9 
        EXPECT_EQ(densityMatrixDiags.size(), 9);
        // We have 2 qubit systems
        const auto finalPopulations = qubitReg->getInformation("<O>").as<std::vector<double>>();
        EXPECT_EQ(finalPopulations.size(), 2);
        // Qubit 1 is in ~|0> state @ t = 1700
        EXPECT_NEAR(finalPopulations[1], 0.0, 0.1);
        // Qubit 0 (control) is near 1.0 
        // (note: there are some oscillation in the control qubit as can be seen in the notebook)    
        EXPECT_NEAR(finalPopulations[0], 1.0, 0.25);
    }
    // ===================================================================================================
}

int main(int argc, char **argv) 
{
  xacc::Initialize(argc, argv);
  ::testing::InitGoogleTest(&argc, argv);
  auto ret = RUN_ALL_TESTS();
  xacc::Finalize();
  return ret;
}