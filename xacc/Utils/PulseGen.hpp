#pragma once

#include <vector>
#include <complex>

namespace QuaC {
    std::vector<std::complex<double>> SquarePulse(size_t in_nbSamples, double in_amplitude = 1.0);
    std::vector<std::complex<double>> GaussianPulse(size_t in_nbSamples, double in_sigma, double in_dt = 1.0, double in_amplitude = 1.0);
    std::vector<std::complex<double>> PulseFunc(const std::string& in_functionString, size_t in_nbSamples, double in_dt = 1.0);
}