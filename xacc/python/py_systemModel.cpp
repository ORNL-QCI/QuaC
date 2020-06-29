#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>
#include "PulseChannelController.hpp"
#include "PulseGen.hpp"
#include "PulseSystemModel.hpp"
#include "xacc_service.hpp"
#include "Pulse.hpp"

namespace py = pybind11;

PYBIND11_MODULE(_pyquaC, m) 
{
    m.doc() = "XACC QuaC pluggin Python wrapper"; 
    m.def(
        "createPulseModel",
        []() -> std::shared_ptr<QuaC::PulseSystemModel> {
            // Get a reference to the global system model service
            auto serviceRef = xacc::getService<QuaC::PulseSystemModel>("global");
            serviceRef->reset();
            return serviceRef;
        },
    "");
    m.def(
        "addPulse",
        [](const std::string& in_pulseName, const std::vector<std::complex<double>>& in_pulseData) {
            auto pulse = std::make_shared<xacc::quantum::Pulse>(in_pulseName);
            std::vector<std::vector<double>> samples;
            samples.reserve(in_pulseData.size());
            for (const auto& dataPoint : in_pulseData)
            {
                samples.emplace_back(std::vector<double>{ dataPoint.real(), dataPoint.imag()});
            }
            pulse->setSamples(samples);
            xacc::contributeService(in_pulseName, pulse);
        },
    "")
    .def("SquarePulse", [](size_t nSamples) {return QuaC::SquarePulse(nSamples);}, "")
    .def("PulseFunc", [](const std::string& in_functionString, size_t in_nbSamples){
        return QuaC::PulseFunc(in_functionString, in_nbSamples);
    })
    .def("PulseFunc", [](const std::string& in_functionString, size_t in_nbSamples, double dt){
        return QuaC::PulseFunc(in_functionString, in_nbSamples, dt);
    })
    .def("GaussianPulse", [](size_t in_nbSamples, double in_sigma) {
        return QuaC::GaussianPulse(in_nbSamples, in_sigma);
    })
    .def("GaussianSquare", [](int in_duration, std::complex<double> in_amp, int in_sigma, int in_width) {
        return QuaC::GaussianSquare(in_duration, in_amp, in_sigma, in_width);
    })
    .def("DragPulse", [](int in_duration, std::complex<double> in_amp, int in_sigma, double in_beta) {
        return QuaC::Drag(in_duration, in_amp, in_sigma, in_beta);
    })
    .def("SlepianPulse", [](size_t in_nbSamples, double in_bW, int in_k) {
        return QuaC::SlepianPulse(in_nbSamples, in_bW, in_k);
    })
    .def("createPulse", [](const std::string& name, const std::string& channel) -> std::shared_ptr<xacc::Instruction> {
        return std::make_shared<xacc::quantum::Pulse>(name, channel);
    });

    py::class_<BackendChannelConfigs>(m, "BackendChannelConfigs")
        .def(py::init<>())
        .def_readwrite("loFregs_dChannels", &BackendChannelConfigs::loFregs_dChannels)
        .def_readwrite("loFregs_uChannels", &BackendChannelConfigs::loFregs_uChannels)
        .def_readwrite("dt", &BackendChannelConfigs::dt)
        .def("addOrReplacePulse", &BackendChannelConfigs::addOrReplacePulse);

    py::class_<QuaC::PulseSystemModel, std::shared_ptr<QuaC::PulseSystemModel>>(m, "PulseSystemModel")
        .def(py::init<>())
        .def("name", (const std::string(QuaC::PulseSystemModel::*)()) & QuaC::PulseSystemModel::name, "Get model name")
        .def("setQubitT1", (void(QuaC::PulseSystemModel::*)(size_t, double)) & QuaC::PulseSystemModel::setQubitT1, "Set T1 of a qubit")
        .def("getQubitT1", (double(QuaC::PulseSystemModel::*)(size_t)) & QuaC::PulseSystemModel::getQubitT1, "Get T1 of a qubit")
        .def("loadHamiltonianJson", (bool(QuaC::PulseSystemModel::*)(const std::string&)) & QuaC::PulseSystemModel::loadHamiltonianJson, "Load Hamiltonian from JSON")
        .def("setChannelConfigs", (void(QuaC::PulseSystemModel::*)(const BackendChannelConfigs&)) & QuaC::PulseSystemModel::setChannelConfigs, "Set backend channel configurations")
        .def("setQubitInitialPopulation", (void(QuaC::PulseSystemModel::*)(size_t, double)) & QuaC::PulseSystemModel::setQubitInitialPopulation, "Set initial state of qubit");
    
}