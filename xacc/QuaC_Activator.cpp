#include "cppmicroservices/BundleActivator.h"
#include "cppmicroservices/BundleContext.h"
#include "cppmicroservices/ServiceProperties.h"
#include "QuaC_Accelerator.hpp"
#include <iostream>
#include "Hamiltonian.hpp"
#include "QuaC_TearDown.hpp"
#include "PulseSystemModel.hpp"
#include "BackendRegistry.hpp"

using namespace cppmicroservices;

class US_ABI_LOCAL QuaC_Activator : public BundleActivator {
public:
  QuaC_Activator() {}

  void Start(BundleContext context) {
    auto acc = std::make_shared<QuaC::QuaC_Accelerator>();
    // Register the QuaC accelerator with the service registry
    context.RegisterService<xacc::Accelerator>(acc);
    context.RegisterService<QuaC::HamiltonianParsingUtil>(std::make_shared<QuaC::HamiltonianParsingUtil>());
    context.RegisterService<xacc::TearDown>(std::make_shared<QuacTearDown>());
    context.RegisterService<QuaC::PulseSystemModel>(std::make_shared<QuaC::PulseSystemModel>("global"));
    context.RegisterService<QuaC::BackendRegistry>(std::make_shared<QuaC::BackendRegistry>());
  }

  void Stop(BundleContext context) {}
};

CPPMICROSERVICES_EXPORT_BUNDLE_ACTIVATOR(QuaC_Activator)
