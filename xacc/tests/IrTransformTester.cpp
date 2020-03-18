#include <gtest/gtest.h>
#include "xacc.hpp"
#include "PulseSystemModel.hpp"

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
    // LO freq = Cavity Freq = 0.0
    channelConfigs.loFregs_dChannels.emplace_back(0.0);      
    systemModel->setChannelConfigs(channelConfigs);
    // Get QuaC accelerator
    auto quaC = xacc::getAccelerator("QuaC", { std::make_pair("system-model", systemModel) });    

    auto xasmCompiler = xacc::getCompiler("xasm");
    auto ir = xasmCompiler->compile(R"(__qpu__ void test(qbit q) {
      X(q[0]);
    })", quaC);

    auto program = ir->getComposite("test");
    auto opt = xacc::getIRTransformation("quantum-control");
    opt->apply(program, quaC);
}


int main(int argc, char **argv) 
{
  xacc::Initialize(argc, argv);
  ::testing::InitGoogleTest(&argc, argv);
  auto ret = RUN_ALL_TESTS();
  xacc::Finalize();
  return ret;
}