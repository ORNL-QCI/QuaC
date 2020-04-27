#include "PulseGen.hpp"
#include <math.h> 
#include "exprtk/exprtk.hpp"
#include <cassert>

using symbol_table_t = exprtk::symbol_table<double>;
using expression_t = exprtk::expression<double>;
using parser_t = exprtk::parser<double>;

namespace QuaC {
    std::vector<std::complex<double>> SquarePulse(size_t in_nbSamples, double in_amplitude)
    {
        const std::vector<std::complex<double>> result(in_nbSamples, in_amplitude);
        return result;
    }

    std::vector<std::complex<double>> GaussianPulse(size_t in_nbSamples, double in_sigma, double in_dt, double in_amplitude)
    {
        std::vector<std::complex<double>> result;
        result.reserve(in_nbSamples);
        for (size_t i = 0; i < in_nbSamples; ++i)
        {
            result.emplace_back(in_amplitude*std::exp(-std::pow(1.0*i*in_dt, 2) / 2.0 / std::pow(in_sigma, 2.0)));
        }
        return result;
    }

    std::vector<std::complex<double>> PulseFunc(const std::string& in_functionString, size_t in_nbSamples, double in_dt)
    {
        std::vector<std::complex<double>> result;
        result.reserve(in_nbSamples);
        expression_t expression;
        parser_t parser;
        symbol_table_t symbol_table;
        double g_time = 0.0;
        symbol_table.add_variable("t", g_time);
        expression.register_symbol_table(symbol_table);
        parser.compile(in_functionString, expression);
        
        for (size_t i = 0; i < in_nbSamples; ++i)
        {
            g_time = i * in_dt;
            result.emplace_back(expression.value());
        }
        return result;
    }

    // Analytical form: 
    // (1) Rising edge:
    // risefall = (duration - width) / 2
    // 0 <= x < risefall
    // f(x) = amp * exp( -(1/2) * (x - risefall/2)^2 / sigma^2) )
    // (2) Square body:
    // risefall <= x < risefall + width
    // f(x) = amp
    // (3) Falling edge:
    // risefall + width <= x < duration
    // f(x) = amp * exp( -(1/2) * (x - (risefall + width)/2)^2 / sigma^2) )
    std::vector<std::complex<double>> GaussianSquare(int in_duration, std::complex<double> in_amp, int in_sigma, int in_width)
    {
        const auto gaussCalc = [&](double in_time, double in_center) {

            const auto x = (in_time-in_center)/(1.0*in_sigma);
            const auto gauss = in_amp*std::exp(-x*x/2.0);
            return gauss;
        };
        
        assert(in_duration > in_width);
        std::vector<std::complex<double>> result;
        result.reserve(in_duration);
        const int risefall = (in_duration - in_width) / 2;
        for (int i = 0; i < risefall; ++i)
        {
            const double center = risefall;
            result.emplace_back(gaussCalc(i, center));
        }
        for (int i = 0; i < in_width; ++i)
        {
            result.emplace_back(in_amp);
        }

        const auto startFall = result.size();
        for (int i = startFall; i < in_duration; ++i)
        {
            result.emplace_back(gaussCalc(i, startFall));
        }
        return result;
    }

    std::vector<std::complex<double>> Drag(int in_duration, std::complex<double> in_amp, int in_sigma, double in_beta)
    {
        std::vector<std::complex<double>> gaussian;
        gaussian.reserve(in_duration);
        const int mu = in_duration/2;
        for (int i = 0; i < in_duration; ++i)
        {
            gaussian.emplace_back(in_amp*std::exp(-std::pow(1.0*i - mu, 2) / 2.0 / std::pow(in_sigma, 2.0)));
        }

        for (int i = 0; i < in_duration; ++i)
        {
            const auto x = gaussian[i];
            static const std::complex<double> I {0.0, 1.0};
            const std::complex<double> deriv_x = I*in_beta*(-(1.0*i - in_duration/2)/(in_sigma*in_sigma))*x;
            gaussian[i] += deriv_x;
        }

        return gaussian;
    }

}