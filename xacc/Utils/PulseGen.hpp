#pragma once

#include <vector>
#include <complex>
#include<armadillo>


namespace QuaC {
    std::vector<std::complex<double>> SquarePulse(size_t in_nbSamples, double in_amplitude = 1.0);
    std::vector<std::complex<double>> GaussianPulse(size_t in_nbSamples, double in_sigma, double in_dt = 1.0, double in_amplitude = 1.0);
    std::vector<std::complex<double>> PulseFunc(const std::string& in_functionString, size_t in_nbSamples, double in_dt = 1.0);
    
    // A square pulse with a Gaussian shaped risefall on either side:
    // PARAMS:
    //  - in_duration (int) – Pulse length in terms of the the sampling period dt.
    //  - in_amp (complex) – The amplitude of the Gaussian and of the square pulse.
    //  - in_sigma (int) – The rise/fall width
    //  - in_width (int) – The duration of the embedded square pulse.
    std::vector<std::complex<double>> GaussianSquare(int in_duration, std::complex<double> in_amp, int in_sigma, int in_width);

    // DRAG pulse:
    std::vector<std::complex<double>> Drag(int in_duration, std::complex<double> in_amp, int in_sigma, double in_beta);

    // Slepian pulse:
    std::vector<double> SlepianPulse(std::vector<double> alpha_vector, size_t in_nbSamples, double in_bW, int in_K);

    // Hanning pulse:
    std::vector<double> HanningPulse(std::vector<double> alpha_vector, size_t in_nbSamples, size_t in_K, double in_T, double in_dt);
}