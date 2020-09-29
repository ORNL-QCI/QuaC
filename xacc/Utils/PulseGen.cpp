#include "PulseGen.hpp"
#include <math.h> 
#include "exprtk/exprtk.hpp"
#include <cassert>
#include <iostream>

using symbol_table_t = exprtk::symbol_table<double>;
using expression_t = exprtk::expression<double>;
using parser_t = exprtk::parser<double>;

namespace {
    template<typename T>
    std::vector<T> arange(T start, T stop, T step = 1) 
    {
        std::vector<T> values;
        for (T value = start; value < stop; value += step)
        {
            values.emplace_back(value);
        }
        return values;
    }

    std::vector<double> genTimePoints(int in_duration, QuaC::SamplerType in_samplerType) 
    {
        switch (in_samplerType)
        {
            case QuaC::SamplerType::Left: return arange<double>(0.0, in_duration);
            case QuaC::SamplerType::Midpoint: return arange(0.5, 0.5 + in_duration);
            case QuaC::SamplerType::Right: return arange<double>(1, in_duration + 1);
        }
        return {};
    }

    std::complex<double> gaussianPoint(double in_time, std::complex<double> in_amp, int in_center, double in_sigma)
    {
        return (in_amp*std::exp(-std::pow(in_time - in_center, 2) / 2.0 / std::pow(in_sigma, 2.0)));
    }
    
    std::vector<std::complex<double>> fixGaussianWidth(const std::vector<std::complex<double>>& in_gauss, const std::complex<double>& in_amp, int in_center, double in_sigma, std::optional<int> in_zeroWidth, bool in_rescaleAmp)
    {
        const int zeroedWidth = in_zeroWidth.value_or(2 * (in_center + 1));
        const auto zeroOffset = gaussianPoint(zeroedWidth / 2, in_amp, 0, in_sigma);
        auto gaussianSamples = in_gauss;
        for (size_t i = 0; i < gaussianSamples.size(); ++i)
        {
            gaussianSamples[i] -= zeroOffset;
        }
        std::complex<double> ampScaleFactor = 1.0;
        if (in_rescaleAmp)
        {
            if (std::abs(in_amp - zeroOffset) != 0.0) 
            {
                ampScaleFactor = in_amp/(in_amp - zeroOffset);
            }
            for (size_t i = 0; i < gaussianSamples.size(); ++i)
            {
                gaussianSamples[i] *= ampScaleFactor;
            }
        }

        return gaussianSamples;
    }
}

namespace QuaC {
    std::vector<std::complex<double>> SquarePulse(size_t in_nbSamples, const std::complex<double>& in_amplitude)
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

    std::vector<std::complex<double>> Drag(int in_duration, std::complex<double> in_amp, int in_sigma, double in_beta, bool in_zeroEnded, SamplerType in_sampler)
    {
        const auto timePoints = genTimePoints(in_duration, in_sampler);
        assert(timePoints.size() == in_duration);
        std::vector<std::complex<double>> gaussian;
        gaussian.reserve(in_duration);
        const int mu = in_duration/2;
        for (const auto& sampleTime :  timePoints)
        {
            const auto x = (sampleTime - mu)/ in_sigma;
            gaussian.emplace_back(in_amp*std::exp(-(x*x/2.0)));
        }

        const std::optional<int> zeroedWidth = in_zeroEnded ? (in_duration + 2) : std::optional<int>();
        const bool rescaleAmp = in_zeroEnded;
        // Fix Gaussian width
        gaussian =  fixGaussianWidth(gaussian, in_amp, mu, in_sigma, zeroedWidth, rescaleAmp);
        
        for (int i = 0; i < timePoints.size(); ++i)
        {
            const auto x = gaussian[i];
            static const std::complex<double> I {0.0, 1.0};
            const std::complex<double> deriv_x = I*in_beta*(-(1.0*timePoints[i] - in_duration/2)/(in_sigma*in_sigma))*x;
            gaussian[i] += deriv_x;
        }

        return gaussian;
    }

    
    // Pseudocode:
    // 1. Create the tridiagonal matrix of size (nbSamples, nbSamples)
    // 2. Solve for the eigenvectors of the tridiagonal system and flip
    //    to be ordered from greatest eigenvalue to least
    // 3. Even orders (k=0,2,4,...) need to have a positive average and
    //    odd orders (k=1,2,3,...) need to start with a positive lobe. 
    // 4. Multiply alpha_vector (size: (in_K, 1)) by the post-processed
    //    eigenVectors matrix (size: (in_nbSamples, in_K)) and sum along rows 
    //    to get final pulse vector   
    std::vector<double> SlepianPulse(std::vector<double> alpha_vector, size_t in_nbSamples, double in_bW, int in_K) 
    {
        // This implementation follows the Discrete Prolate Spheroidal Sequences 
        // algorithm as outlined in: Percival DB, Walden WT. Spectral Analysis for 
        // Physical Applications: Multitaper and Conventional Univariate Techniques.
        // Cambridge University Press; 1993.

        if (in_nbSamples % 2 != 0)
        {
            std::cout << " Requested odd number of samples : " << in_nbSamples << std::endl;
            std::cout << " Using " << in_nbSamples + 1 << " samples instead " << std::endl;
            size_t in_nbSamples = in_nbSamples + 1;
        }
        
        // Step 1
        arma::vec alphas(alpha_vector);
        arma::SpMat<double> tridiag = arma::zeros<arma::SpMat<double>>(in_nbSamples, in_nbSamples);
        for (size_t i = 0; i < in_nbSamples; i++) 
        {
            for (size_t j = 0; j < in_nbSamples; j++) 
            {
            if (j == i - 1) 
            {
                tridiag(i, j) = 0.5 * i * (in_nbSamples - i);
            }
            if (j == i) 
            {
                tridiag(i, j) = ((0.5 * (in_nbSamples - 1)) - i) *
                                ((0.5 * (in_nbSamples - 1)) - i) *
                                cos(2 * M_PI * in_bW);
            }
            if (j == i + 1) 
            {
                tridiag(i, j) = 0.5 * (i + 1) * (in_nbSamples - 1 - i);
            }
            }
        }

        // Step 2
        arma::vec eigenValues;
        arma::mat eigenVectors_raw;
        arma::eigs_sym(eigenValues, eigenVectors_raw, tridiag, in_K);
        arma::mat eigenVectors = arma::reverse(eigenVectors_raw);

        // Step 3
        auto thresh = std::max(1e-7, 1. / in_nbSamples);
        for (size_t i = 0; i < eigenVectors.n_cols; i++)
        {
            // If it's an even order (k=0,2,4,..), check if it has a positive average.
            // If not, multiply by -1. Per
            if (i % 2 == 0)
            {
                double sum = arma::accu(eigenVectors.col(i));
                if (sum < 0.)
                {
                    for (size_t j = 0; j < eigenVectors.n_rows; ++j)
                    {
                    eigenVectors(j, i) =  eigenVectors(j, i) * -1.;
                    }
                }
            } else {
                // For anti-symmetric orders (k=1,3,5,...) need to begin with a positive lobe
                arma::vec column = eigenVectors.col(i);
                arma::mat ww = column % column;
                arma::vec find_thresh = column.elem( arma::find(ww > thresh) ) ;
                double first_idx = find_thresh.at(0);
                if (first_idx < 0.)
                {
                    for (size_t j = 0; j < eigenVectors.n_rows; ++j)
                    {
                    eigenVectors(j, i) =  eigenVectors(j, i) * -1.;
                    }
                }
            }
        }
        
        // Step 4
        for (size_t i = 0; i < eigenVectors.n_cols; ++i)
        {
          for (size_t j = 0; j < eigenVectors.n_rows; ++j)
          {
            double weight = alphas(i);
            eigenVectors(j, i) = weight * eigenVectors(j, i); 
          }
        } 

        arma::colvec result = arma::sum(eigenVectors, 1);
        return arma::conv_to< std::vector<double> >::from(result);
    }

    // Creating pulse as a weighted expansion in the Hanning basis. 
    // Analytical Form:
    // Omega = \sum{k=1}^{K} \alpha_{k} * (1 - cos(2 * \pi * k * t) / T)
    std::vector<double> HanningPulse(std::vector<double> alpha_vector, size_t in_nbSamples, size_t in_K, size_t in_T)
    {    
        arma::vec order;
        arma::vec result(in_T, arma::fill::zeros);
        for (size_t k = 0; k < in_K; k++) 
        {
            order.zeros(in_T) ;
            for (size_t t = 0; t < in_T; t++) 
            {
                order(t) = alpha_vector.at(k) * (1 - cos((2 * M_PI * k * t) / in_T)) ; 
            }
            result = result + order ;
        }
        return arma::conv_to< std::vector<double> >::from(result);
    }
}