#pragma once

#include <memory>
#include "Identifiable.hpp"

namespace xacc { namespace quantum {
    class Pulse;
}} 

namespace QuaC {
// Convert a parametric pulse (shape + params) into concrete sample pulse instruction for simulation.
class ParametricPulses : public xacc::Identifiable
{
public:
    const std::string name() const override { return "default"; }
    const std::string description() const override { return ""; }
    std::shared_ptr<xacc::quantum::Pulse> generatePulse(const std::string &in_shape, double dt, const std::string &in_paramsJson);
};
}