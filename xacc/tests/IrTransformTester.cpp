#include <gtest/gtest.h>
#include "xacc.hpp"
#include "PulseSystemModel.hpp"
#include "CommonGates.hpp"

TEST(IrTransformTester, checkSimple) 
{
    // Note: this is Hamiltonian is in the rotating frame
    const std::string singleQubitHamiltonianJson = R"#(
    {
        "description": "One-qubit Hamiltonian (used by GOAT)",
        "h_latex": "",
        "h_str": ["omega0*Z0", "omegaa*X0||D0"],
        "osc": {},
        "qub": {
            "0": 2
        },
        "vars": {
            "omega0": 0.0,
            "omegaa": 0.062832
        }
    }
    )#";

    auto systemModel = std::make_shared<QuaC::PulseSystemModel>();
    systemModel->loadHamiltonianJson(singleQubitHamiltonianJson);    
    
    BackendChannelConfigs channelConfigs;
    channelConfigs.dt = 1.0;
    // Rotating frame
    channelConfigs.loFregs_dChannels.emplace_back(0.0);      
    systemModel->setChannelConfigs(channelConfigs);
    // Get QuaC accelerator
    auto quaC = xacc::getAccelerator("QuaC", { 
      std::make_pair("system-model", systemModel), 
      std::make_pair("shots", 10000) 
    });    

    auto xasmCompiler = xacc::getCompiler("xasm");
    // A *complex* way to do X gate (H-Z-H)
    auto ir = xasmCompiler->compile(R"(__qpu__ void test(qbit q) {
      H(q[0]);
      Z(q[0]);
      H(q[0]);
    })", quaC);

    auto program = ir->getComposite("test");

    // Pulse IR transformation configs:
    const std::vector<double> initParams { 8.0 };
    const double tMax = 100;
    xacc::HeterogeneousMap configs {
      std::make_pair("method", "GOAT"),
      std::make_pair("control-params", std::vector<std::string> { "sigma" }),
      // Gaussian pulse
      std::make_pair("control-funcs", std::vector<std::string> { "exp(-t^2/(2*sigma^2))" }),
      // Initial params
      std::make_pair("initial-parameters", initParams),
      std::make_pair("max-time", tMax)
    };

    auto opt = xacc::getIRTransformation("quantum-control");
    if (opt) 
    {
      opt->apply(program, quaC, configs);
      // This should be a pulse now
      EXPECT_EQ(program->nInstructions(), 1);
      EXPECT_TRUE(program->getInstruction(0)->isAnalog());
      // Simulate the optimized pulse program
      auto qubitReg = xacc::qalloc(1);    
      // Add a measurement to get bit count statistics
      auto meas = std::make_shared<xacc::quantum::Measure>(0);
      program->addInstruction(meas);
      quaC->execute(qubitReg, program);
      // Should be an X gate
      qubitReg->print();
      // Check via bitstrings
      const double prob0 = qubitReg->computeMeasurementProbability("0");
      const double prob1 = qubitReg->computeMeasurementProbability("1");
      EXPECT_NEAR(prob0, 0.0, 0.01);
      EXPECT_NEAR(prob1, 1.0, 0.01);
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