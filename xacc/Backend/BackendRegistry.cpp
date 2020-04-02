#include "BackendRegistry.hpp"
#include "PulseSystemModel.hpp"
#include "DefaultBackend2Q.hpp"

namespace {
std::vector<std::shared_ptr<QuaC::PulseSystemModel>>& getSystemModels()
{
    static std::vector<std::shared_ptr<QuaC::PulseSystemModel>> g_allModels {
        std::make_shared<QuaC::Default2Q>()
    };
    return g_allModels;
}
}

namespace QuaC {
std::shared_ptr<PulseSystemModel> BackendRegistry::getSystemModel(const std::string& in_name)
{
    for (const auto& model : getSystemModels())
    {
        if (model->name() == in_name)
        {
            return model;
        }
    }

    return nullptr;
}
}