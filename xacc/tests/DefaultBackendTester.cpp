#include <gtest/gtest.h>
#include "xacc.hpp"

TEST(DefaultBackendTester, checkGateX) 
{
    const int NB_SHOTS = 10000;
    {
        auto qpu = xacc::getAccelerator("QuaC:Default2Q", {std::make_pair("shots", NB_SHOTS)});    
    
        auto xasmCompiler = xacc::getCompiler("xasm");
        // X gate:
        auto ir = xasmCompiler->compile(R"(__qpu__ void test1(qbit q) {
            X(q[0]);
            Measure(q[0]);
        })", qpu);

        auto program = ir->getComposite("test1");
        auto qubitReg = xacc::qalloc(2);    
        qpu->execute(qubitReg, program);
        qubitReg->print();
        const double prob0 = qubitReg->computeMeasurementProbability("0");
        const double prob1 = qubitReg->computeMeasurementProbability("1");
        // Should be in |1> state
        EXPECT_NEAR(prob0, 0.0, 0.1);
        EXPECT_NEAR(prob1, 1.0, 0.1);  
    }
    
    // {
    //     auto qpu = xacc::getAccelerator("QuaC:Default2Q", {std::make_pair("shots", NB_SHOTS)});    
    
    //     auto xasmCompiler = xacc::getCompiler("xasm");
    //     // H gate:
    //     auto ir = xasmCompiler->compile(R"(__qpu__ void test2(qbit q) {
    //         H(q[0]);
    //         Measure(q[0]);
    //     })", qpu);

    //     auto program = ir->getComposite("test2");
    //     auto qubitReg = xacc::qalloc(2);    
    //     qpu->execute(qubitReg, program);
    //     const double prob0 = qubitReg->computeMeasurementProbability("0");
    //     const double prob1 = qubitReg->computeMeasurementProbability("1");
    //     // Should be in 50-50 state
    //     EXPECT_NEAR(prob0, 0.5, 0.1);
    //     EXPECT_NEAR(prob1, 0.5, 0.1);  
    // }
    // {
    //     auto qpu = xacc::getAccelerator("QuaC:Default2Q", {std::make_pair("shots", NB_SHOTS)});    
    
    //     auto xasmCompiler = xacc::getCompiler("xasm");
    //     // X gate:
    //     auto ir = xasmCompiler->compile(R"(__qpu__ void test3(qbit q) {
    //         X(q[1]);
    //         Measure(q[1]);
    //     })", qpu);

    //     auto program = ir->getComposite("test3");
    //     auto qubitReg = xacc::qalloc(2);    
    //     qpu->execute(qubitReg, program);
    //     const double prob0 = qubitReg->computeMeasurementProbability("0");
    //     const double prob1 = qubitReg->computeMeasurementProbability("1");
    //     // Should be in |1> state
    //     EXPECT_NEAR(prob0, 0.0, 0.1);
    //     EXPECT_NEAR(prob1, 1.0, 0.1);  
    // }
    
    // {
    //     auto qpu = xacc::getAccelerator("QuaC:Default2Q", {std::make_pair("shots", NB_SHOTS)});    
    
    //     auto xasmCompiler = xacc::getCompiler("xasm");
    //     // H gate:
    //     auto ir = xasmCompiler->compile(R"(__qpu__ void test4(qbit q) {
    //         H(q[1]);
    //         Measure(q[1]);
    //     })", qpu);

    //     auto program = ir->getComposite("test4");
    //     auto qubitReg = xacc::qalloc(2);    
    //     qpu->execute(qubitReg, program);
    //     const double prob0 = qubitReg->computeMeasurementProbability("0");
    //     const double prob1 = qubitReg->computeMeasurementProbability("1");
    //     // Should be in 50-50 state
    //     EXPECT_NEAR(prob0, 0.5, 0.1);
    //     EXPECT_NEAR(prob1, 0.5, 0.1);  
    // }
}

int main(int argc, char **argv) 
{
  xacc::Initialize(argc, argv);
  ::testing::InitGoogleTest(&argc, argv);
  auto ret = RUN_ALL_TESTS();
  xacc::Finalize();
  return ret;
}
