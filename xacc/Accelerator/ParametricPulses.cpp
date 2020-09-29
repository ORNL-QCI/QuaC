#include "ParametricPulses.hpp"
#include "PulseGen.hpp"
#include "json.hpp"
#include "xacc.hpp"
#include "Pulse.hpp"

namespace QuaC {
std::shared_ptr<xacc::quantum::Pulse> ParametricPulses::generatePulse(const std::string &in_shape, const std::string &in_paramsJson) 
{
    static const std::vector<std::string> SUPPORTED_PULSE_SHAPES { "gaussian", "gaussian_square", "drag", "constant" };
    if (!xacc::container::contains(SUPPORTED_PULSE_SHAPES, in_shape)) 
    {
        xacc::error("Pulse shape '" + in_shape + "' is not supported.");
        return nullptr;
    }
    
    const auto formatPulseSampleVec = [](const std::vector<std::complex<double>>& in_pulseData) -> std::vector<std::vector<double>> {
        std::vector<std::vector<double>> samples;
        samples.reserve(in_pulseData.size());
        for (const auto& dataPoint : in_pulseData)
        {
            samples.emplace_back(std::vector<double>{ dataPoint.real(), dataPoint.imag()});
        }
        return samples;
    };

    if (in_shape == "gaussian") 
    {
        static int pulseIdCounter = 0;
        auto pulseInst = std::make_shared<xacc::quantum::Pulse>("gaussian_" + std::to_string(pulseIdCounter++));
        auto j = nlohmann::json::parse(in_paramsJson);   
        auto ampAsVec = j["amp"].get<std::vector<double>>();
        assert(ampAsVec.size() == 2);
        auto duration = j["duration"].get<int>();
        auto sigma = j["sigma"].get<int>();
        const std::complex<double> amp { ampAsVec[0], ampAsVec[1]};
        // Beta = 0.0 (no derivative part)
        const auto pulseSamples = QuaC::Drag(duration, amp, sigma, 0.0);
        pulseInst->setSamples(formatPulseSampleVec(pulseSamples));
        return pulseInst;
    }
    
    if (in_shape == "gaussian_square") 
    {
        static int pulseIdCounter = 0;
        auto pulseInst = std::make_shared<xacc::quantum::Pulse>("gaussian_square_" + std::to_string(pulseIdCounter++));
        auto j = nlohmann::json::parse(in_paramsJson);   
        auto ampAsVec = j["amp"].get<std::vector<double>>();
        assert(ampAsVec.size() == 2);
        auto width = j["width"].get<int>();
        auto duration = j["duration"].get<int>();
        assert(width < duration);
        auto sigma = j["sigma"].get<int>();
        const std::complex<double> amp { ampAsVec[0], ampAsVec[1]};
        const auto pulseSamples = QuaC::GaussianSquare(duration, amp, sigma, width);
        pulseInst->setSamples(formatPulseSampleVec(pulseSamples));
        return pulseInst;
    }

    if (in_shape == "drag") 
    {
        static int pulseIdCounter = 0;
        auto pulseInst = std::make_shared<xacc::quantum::Pulse>("drag_" + std::to_string(pulseIdCounter++));
        auto j = nlohmann::json::parse(in_paramsJson);   
        auto ampAsVec = j["amp"].get<std::vector<double>>();
        assert(ampAsVec.size() == 2);
        auto beta = j["beta"].get<double>();
        auto duration = j["duration"].get<int>();
        auto sigma = j["sigma"].get<int>();
        const std::complex<double> amp { ampAsVec[0], ampAsVec[1]};
        const auto pulseSamples = QuaC::Drag(duration, amp, sigma, beta, true, QuaC::SamplerType::Midpoint);
        pulseInst->setSamples(formatPulseSampleVec(pulseSamples));
        return pulseInst;
    }

    if (in_shape == "constant") 
    {
        static int pulseIdCounter = 0;
        auto pulseInst = std::make_shared<xacc::quantum::Pulse>("constant_" + std::to_string(pulseIdCounter++));
        auto j = nlohmann::json::parse(in_paramsJson);   
        auto ampAsVec = j["amp"].get<std::vector<double>>();
        assert(ampAsVec.size() == 2);
        auto duration = j["duration"].get<int>();
        const std::complex<double> amp { ampAsVec[0], ampAsVec[1]};
        const auto pulseSamples = QuaC::SquarePulse(duration, amp);
        pulseInst->setSamples(formatPulseSampleVec(pulseSamples));
        return pulseInst;
    }

    return nullptr;
}
}