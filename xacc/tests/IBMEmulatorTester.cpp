#include <gtest/gtest.h>
#include "xacc.hpp"
#include "xacc_service.hpp"


TEST(IBMEmulatorTester, checkSimple) 
{
    xacc::set_verbose(true);
    auto qpu = xacc::getAccelerator("QuaC:ibmq_armonk"); 
    auto buffer = xacc::qalloc(1);
    auto provider = xacc::getService<xacc::IRProvider>("quantum");
    auto f = provider->createComposite("tmp");
    auto m = provider->createInstruction("Measure", {0});
    auto h = provider->createInstruction("H", {0});
    f->addInstruction(h);
    f->addInstruction(m);
    qpu->execute(buffer, f);
    buffer->print();
}

int main(int argc, char **argv) 
{
  xacc::Initialize(argc, argv);
  ::testing::InitGoogleTest(&argc, argv);
  auto ret = RUN_ALL_TESTS();
  xacc::Finalize();
  return ret;
}
