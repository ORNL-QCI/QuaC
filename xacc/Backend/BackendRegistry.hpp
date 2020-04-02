#pragma once
#include "xacc.hpp"

namespace QuaC {
class PulseSystemModel;
class BackendRegistry : public xacc::Identifiable
{
public:
    virtual const std::string name() const override { return "default"; }
    virtual const std::string description() const { return ""; }
    // Null if name is invalid.
    std::shared_ptr<PulseSystemModel> getSystemModel(const std::string& in_name);
};
}