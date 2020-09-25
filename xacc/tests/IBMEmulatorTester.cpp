#include <gtest/gtest.h>
#include "xacc.hpp"


TEST(IBMEmulatorTester, checkSimple) 
{
    auto qpu = xacc::getAccelerator("QuaC:ibmq_bogota"); 
}

int main(int argc, char **argv) 
{
  xacc::Initialize(argc, argv);
  ::testing::InitGoogleTest(&argc, argv);
  auto ret = RUN_ALL_TESTS();
  xacc::Finalize();
  return ret;
}
