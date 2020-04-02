#include <gtest/gtest.h>
#include "xacc.hpp"


TEST(DefaultBackendTester, checkSimple) 
{
    auto qpu = xacc::getAccelerator("QuaC:Default2Q");    

}

int main(int argc, char **argv) 
{
  xacc::Initialize(argc, argv);
  ::testing::InitGoogleTest(&argc, argv);
  auto ret = RUN_ALL_TESTS();
  xacc::Finalize();
  return ret;
}
