#include <gtest/gtest.h>
#include "xacc.hpp"

namespace {
    // Tolerance epsilon for numerical errors (5%)
    const double TOL_EPS = 0.05;
}

TEST(DefaultBackendTester, checkSingleGateQ0) 
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
        EXPECT_NEAR(prob0, 0.0, TOL_EPS);
        EXPECT_NEAR(prob1, 1.0, TOL_EPS);  
    }
    
    {
        auto qpu = xacc::getAccelerator("QuaC:Default2Q", {std::make_pair("shots", NB_SHOTS)});    
    
        auto xasmCompiler = xacc::getCompiler("xasm");
        // H gate:
        auto ir = xasmCompiler->compile(R"(__qpu__ void test2(qbit q) {
            H(q[0]);
            Measure(q[0]);
        })", qpu);

        auto program = ir->getComposite("test2");
        auto qubitReg = xacc::qalloc(2);    
        qpu->execute(qubitReg, program);
        qubitReg->print();

        const double prob0 = qubitReg->computeMeasurementProbability("0");
        const double prob1 = qubitReg->computeMeasurementProbability("1");
        // Should be in 50-50 state
        EXPECT_NEAR(prob0, 0.5, TOL_EPS);
        EXPECT_NEAR(prob1, 0.5, TOL_EPS);  
    }

    {
        auto qpu = xacc::getAccelerator("QuaC:Default2Q", {std::make_pair("shots", NB_SHOTS)});    
        // Test Hadamard-Hadamard is equivalent to Identity
        auto xasmCompiler = xacc::getCompiler("xasm");
        // H gate:
        auto ir = xasmCompiler->compile(R"(__qpu__ void test3(qbit q) {
            H(q[0]);
            H(q[0]);
            Measure(q[0]);
        })", qpu);

        auto program = ir->getComposite("test3");
        auto qubitReg = xacc::qalloc(2);    
        qpu->execute(qubitReg, program);
        qubitReg->print();

        const double prob0 = qubitReg->computeMeasurementProbability("0");
        const double prob1 = qubitReg->computeMeasurementProbability("1");
        EXPECT_NEAR(prob0, 1.0, TOL_EPS);
        EXPECT_NEAR(prob1, 0.0, TOL_EPS);  
    }

    {
        auto qpu = xacc::getAccelerator("QuaC:Default2Q", {std::make_pair("shots", NB_SHOTS)});    
        auto xasmCompiler = xacc::getCompiler("xasm");
        auto ir = xasmCompiler->compile(R"(__qpu__ void testRx(qbit q, double t) {
            Rx(q[0], t);
            Measure(q[0]);
        })", qpu);

        auto program = ir->getComposite("testRx");

        const auto angles = xacc::linspace(-xacc::constants::pi, xacc::constants::pi, 20);
        for (const auto& a : angles) 
        {
            auto buffer = xacc::qalloc(2);
            auto evaled = program->operator()({a});
            qpu->execute(buffer, evaled);

            const double prob0 = buffer->computeMeasurementProbability("0");
            const double expectedResult = std::cos(a/2.0)*std::cos(a/2.0);
            std::cout << "Angle = " << a << ": Prob(0) = " << prob0 << " vs. Expected = " << expectedResult << "\n";
            EXPECT_NEAR(prob0, expectedResult, TOL_EPS);
        }
    }
}

TEST(DefaultBackendTester, checkSingleGateQ1) 
{
    const int NB_SHOTS = 10000;
    {
        auto qpu = xacc::getAccelerator("QuaC:Default2Q", {std::make_pair("shots", NB_SHOTS)});    
    
        auto xasmCompiler = xacc::getCompiler("xasm");
        // X gate:
        auto ir = xasmCompiler->compile(R"(__qpu__ void test4(qbit q) {
            X(q[1]);
            Measure(q[1]);
        })", qpu);

        auto program = ir->getComposite("test4");
        auto qubitReg = xacc::qalloc(2);    
        qpu->execute(qubitReg, program);
        qubitReg->print();
        const double prob0 = qubitReg->computeMeasurementProbability("0");
        const double prob1 = qubitReg->computeMeasurementProbability("1");
        // Should be in |1> state
        EXPECT_NEAR(prob0, 0.0, TOL_EPS);
        EXPECT_NEAR(prob1, 1.0, TOL_EPS);  
    }
    
    {
        auto qpu = xacc::getAccelerator("QuaC:Default2Q", {std::make_pair("shots", NB_SHOTS)});    
    
        auto xasmCompiler = xacc::getCompiler("xasm");
        // H gate:
        auto ir = xasmCompiler->compile(R"(__qpu__ void test5(qbit q) {
            H(q[1]);
            Measure(q[1]);
        })", qpu);

        auto program = ir->getComposite("test5");
        auto qubitReg = xacc::qalloc(2);    
        qpu->execute(qubitReg, program);
        qubitReg->print();

        const double prob0 = qubitReg->computeMeasurementProbability("0");
        const double prob1 = qubitReg->computeMeasurementProbability("1");
        // Should be in 50-50 state
        EXPECT_NEAR(prob0, 0.5, TOL_EPS);
        EXPECT_NEAR(prob1, 0.5, TOL_EPS);  
    }

    {
        auto qpu = xacc::getAccelerator("QuaC:Default2Q", {std::make_pair("shots", NB_SHOTS)});    
        // Test Hadamard-Hadamard is equivalent to Identity
        auto xasmCompiler = xacc::getCompiler("xasm");
        // H gate:
        auto ir = xasmCompiler->compile(R"(__qpu__ void test6(qbit q) {
            H(q[1]);
            H(q[1]);
            Measure(q[0]);
        })", qpu);

        auto program = ir->getComposite("test6");
        auto qubitReg = xacc::qalloc(2);    
        qpu->execute(qubitReg, program);
        qubitReg->print();

        const double prob0 = qubitReg->computeMeasurementProbability("0");
        const double prob1 = qubitReg->computeMeasurementProbability("1");
        EXPECT_NEAR(prob0, 1.0, TOL_EPS);
        EXPECT_NEAR(prob1, 0.0, TOL_EPS);  
    }

    {
        auto qpu = xacc::getAccelerator("QuaC:Default2Q", {std::make_pair("shots", NB_SHOTS)});    
        auto xasmCompiler = xacc::getCompiler("xasm");
        auto ir = xasmCompiler->compile(R"(__qpu__ void testRx1(qbit q, double t) {
            Rx(q[1], t);
            Measure(q[1]);
        })", qpu);

        auto program = ir->getComposite("testRx1");

        const auto angles = xacc::linspace(-xacc::constants::pi, xacc::constants::pi, 20);
        for (const auto& a : angles) 
        {
            auto buffer = xacc::qalloc(2);
            auto evaled = program->operator()({a});
            qpu->execute(buffer, evaled);

            const double prob0 = buffer->computeMeasurementProbability("0");
            const double expectedResult = std::cos(a/2.0)*std::cos(a/2.0);
            std::cout << "Angle = " << a << ": Prob(0) = " << prob0 << " vs. Expected = " << expectedResult << "\n";
            EXPECT_NEAR(prob0, expectedResult, TOL_EPS);
        }
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
