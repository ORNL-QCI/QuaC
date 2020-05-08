#include "Functor.hpp"
#include <chrono>
#include <iomanip>
#include <fstream>
#include <cassert>
#include <unistd.h>
#include <limits.h>

extern "C" {
#include "interface_xacc_ir.h"
}

namespace {
    std::string writeTimesteppingDataToCsv(const std::string& in_fileName, const TSData* const in_tsData, int in_nbSteps, int in_nbQubits)
    {
        if (in_nbSteps < 1)
        {
            return "";
        }

        const auto stringEndsWith = [](const std::string& in_string, const std::string& in_ending) {
            if (in_ending.size() > in_string.size()) 
            {
            return false;
            }

            return std::equal(in_ending.rbegin(), in_ending.rend(), in_string.rbegin());   
        };

        const auto getCurrentTimeString = [](){
            const auto currentTime = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
            std::stringstream ss;
            ss << std::put_time(std::localtime(&currentTime), "%Y%m%d_%X");
            return ss.str();
        };

        std::ofstream outputFile;
        std::string fileName = stringEndsWith(in_fileName, ".csv") ?  in_fileName.substr(0, in_fileName.size() - 4) : in_fileName;
        // Add a timestamp to prevent duplicate
        fileName += "_";
        fileName += getCurrentTimeString();
        fileName += ".csv";
        outputFile.open(fileName, std::ofstream::out);

        if(!outputFile.is_open())
        {
            std::cout << "Cannot open CSV file '" << fileName << "' for writing!\n";
            return "";
        }

        const auto nbChannels = in_tsData[0].nbChannels;      
        const auto nbPopulations = in_tsData[0].nbPops;
        // Header
        outputFile << "Time, "; 
        for(int j = 0; j < nbChannels; ++j)
        {
            outputFile << "Channel[" << j << "], ";
        }     
        for(int j = 0; j < nbPopulations; ++j)
        {
            outputFile << "Population[" << j << "], ";
        }
        
        for(int j = 0; j < in_nbQubits; ++j)
        {
            outputFile << "<X[" << j << "]>, <Y[" << j << "]>, <Z[" << j << "]>, ";
        }

        outputFile << "\n";

        // Data
        for (int i = 0; i < in_nbSteps; ++i)
        {
            const TSData dataAtStep = in_tsData[i];
            // First column is time:
            outputFile << dataAtStep.time << ", ";
            // Channel data columns
            for(int j = 0; j < nbChannels; ++j)
            {
            outputFile << dataAtStep.channelData[j] << ", ";
            }
            // Population columns
            for(int j = 0; j < nbPopulations; ++j)
            {
            outputFile << dataAtStep.populations[j] << ", ";
            }
            
            // Pauli expectations
            for(int j = 0; j < in_nbQubits; ++j)
            {
            for (int ii = 0; ii < 3; ++ii)
            {
                outputFile << dataAtStep.pauliExpectations[3*j + ii] << ", ";
            }
            }
            
            outputFile << "\n";
        }
        outputFile.close();

        char currentDirectory[PATH_MAX];
        if (getcwd(currentDirectory, sizeof(currentDirectory)) != NULL) 
        {
            fileName = std::string(currentDirectory) + "/" + fileName;
            std::cout << "Time-stepping data is written to file '" << fileName << "'\n";
            return fileName;
        } 

        return "";
    } 
}

InitializeFunctor::InitializeFunctor(int in_nbQubit, 
                const std::vector<int>& in_qbitDims, 
                const std::unordered_map<int, double>& in_qbitDecays, 
                const std::unordered_map<int, double>& in_qbitInitialPopulations,
                bool in_verbose,
                double in_monitorDt):
    m_nbQubit(in_nbQubit),
    m_qbitDims(in_qbitDims),
    m_qbitDecays(m_nbQubit, 0.0),
    m_qbitInitialPops(m_nbQubit, 0.0),
    m_verbose(in_verbose),
    m_monitorDt(in_monitorDt)  
{
    for (const auto& kv : in_qbitDecays) 
    {
        m_qbitDecays[kv.first] = kv.second;
    }
    for (const auto& kv : in_qbitInitialPopulations) 
    {
        m_qbitInitialPops[kv.first] = kv.second;
    }
}

void InitializeFunctor::execute(SerializationType* out_result) 
{
    if (m_verbose)
    {
        XACC_QuaC_SetLogVerbosity(DEBUG);
    }

    XACC_QuaC_Initialize(m_nbQubit, m_qbitDims.data());
    for (int i = 0; i < m_nbQubit; ++i)
    {
        XACC_QuaC_AddQubitDecay(i, m_qbitDecays[i]);
    }
    for (int i = 0; i < m_nbQubit; ++i)
    {
        XACC_QuaC_SetInitialPopulation(i, m_qbitInitialPops[i]);
    }    

    if (m_monitorDt > 0.0)
    {
        XACC_QuaC_EnableTimeSteppingMonitor(m_monitorDt);
    } 
}

AddHamiltonianTerm::AddHamiltonianTerm(const std::complex<double>& in_coeff, 
                                    const std::vector<std::pair<std::string, int>>& in_ops, 
                                    int in_channelId):
    m_channelId(in_channelId),
    m_coeff(in_coeff),
    m_ops(in_ops)
    {}

void AddHamiltonianTerm::execute(SerializationType* out_result)
{
    assert(!m_ops.empty() && m_ops.size()<=2);
    if (m_channelId >= 0)
    {
        if (m_ops.size() == 1)
        {
            XACC_QuaC_AddTimeDependentHamiltonianTerm1(
                m_ops[0].first.c_str(), m_ops[0].second,
                m_channelId,
                m_coeff.real());
        }
        else
        {
            XACC_QuaC_AddTimeDependentHamiltonianTerm2(
                m_ops[0].first.c_str(), m_ops[0].second,
                m_ops[1].first.c_str(), m_ops[1].second,
                m_channelId,
                m_coeff.real());
        }
    }
    else
    {
        if (m_ops.size() == 1)
        {
            XACC_QuaC_AddConstHamiltonianTerm1(m_ops[0].first.c_str(), m_ops[0].second, 
                                                { m_coeff.real() , m_coeff.imag()});    
        }
        else
        {
            XACC_QuaC_AddConstHamiltonianTerm2(m_ops[0].first.c_str(), m_ops[0].second, 
                                                m_ops[1].first.c_str(), m_ops[1].second,
                                                { m_coeff.real() , m_coeff.imag()});    
        }
    }
}

AddGateU3::AddGateU3(int in_qubitIdx, double in_theta, double in_phi, double in_lambda, double in_startTime):
    m_qubitIdx(in_qubitIdx),
    m_theta(in_theta),
    m_phi(in_phi),
    m_lambda(in_lambda),
    m_startTime(in_startTime)
{}

void AddGateCNOT::execute(SerializationType* out_result)
{
    XACC_QuaC_AddCnot(m_ctrlIdx, m_targetIdx, m_startTime);
}

AddGateCNOT::AddGateCNOT(int in_ctrlIdx, int in_targetIdx, double in_startTime):
    m_ctrlIdx(in_ctrlIdx),
    m_targetIdx(in_targetIdx),
    m_startTime(in_startTime)
{}

void AddGateU3::execute(SerializationType* out_result)
{
    XACC_QuaC_AddDigitalInstructionU3(m_qubitIdx, m_theta, m_phi, m_lambda, m_startTime);
}

StartTimestepping::StartTimestepping(const PulseChannelController& in_pulseDataProvider, double in_dt, double in_stopTime, int in_stepMax, bool in_adaptive):
    m_pulseDataProvider(in_pulseDataProvider),
    m_dt(in_dt),
    m_stopTime(in_stopTime),
    m_stepMax(in_stepMax),
    m_adaptive(in_adaptive)
{}

void StartTimestepping::execute(SerializationType* out_result)
{
    if (!m_adaptive)
    {
        XACC_QuaC_DisableAdaptiveTimestepping();
    }
    
    double* results = nullptr;
    TSData* tsData;  
    int nbSteps = 0;
    const auto resultSize = XACC_QuaC_RunPulseSim(
                                reinterpret_cast<PulseChannelProvider*>(&m_pulseDataProvider), 
                                m_dt, m_stopTime, m_stepMax, &results, &nbSteps, &tsData);

    // Always retrieve the diagonal elements of the density matrix after finishing the simulation
    ComplexCoefficient* diagElems = nullptr;
    const auto nbDiagonalElems = XACC_QuaC_GetDensityMatrixDiagElements(&diagElems);
    
    // There is no point running time stepping w/o getting the result.
    assert(out_result != nullptr);
    
    SimResult finalResult;
    // If there are time-stepping data recorded
    if (nbSteps > 0)
    {
        finalResult.csvFileName = writeTimesteppingDataToCsv("output", tsData, nbSteps, resultSize);
    }
    
    // Population (occupation expectation) for each qubit
    for (int i = 0; i < resultSize; ++i)
    {
        finalResult.finalPopulations.emplace_back(results[i]);
    }

    // Diagonal Elements:
    finalResult.dmDiagElems.reserve(nbDiagonalElems);
    for (int i = 0; i < nbDiagonalElems; ++i)
    {
        finalResult.dmDiagElems.emplace_back(diagElems[i].real, diagElems[i].imag);
    }

    finalResult.tsData.reserve(nbSteps);
    for (int i = 0; i < nbSteps; ++i)
    {
        const TSData dataAtStep = tsData[i];
        TimeSteppingData data;
        data.time = dataAtStep.time;
        // Channel data columns
        for(int j = 0; j < dataAtStep.nbChannels; ++j)
        {
            data.channelData.emplace_back(dataAtStep.channelData[j]);
        }
        // Population columns
        for(int j = 0; j < dataAtStep.nbPops; ++j)
        {
            data.populations.emplace_back(dataAtStep.populations[j]);
        }
        
        // Pauli expectations
        for(int j = 0; j < resultSize; ++j)
        {
            for (int ii = 0; ii < 3; ++ii)
            {
                data.pauliExpectations.emplace_back(dataAtStep.pauliExpectations[3*j + ii]);
            }
        }

        finalResult.tsData.emplace_back(std::move(data));
    }
    free(results);
    free(diagElems);
    SerializationOutputDataType outArchive(*out_result); 
    outArchive(finalResult); 
}

CalculateBipartiteConcurrence::CalculateBipartiteConcurrence(size_t in_q1, size_t in_q2):
    m_qubit1(in_q1),
    m_qubit2(in_q2)
{}

void CalculateBipartiteConcurrence::execute(SerializationType* out_result)
{
    assert(out_result != nullptr);
    const double concurrentResult = XACC_QuaC_CalcConcurrence(m_qubit1, m_qubit2);
    SerializationOutputDataType outArchive(*out_result); 
    outArchive(concurrentResult); 
}

GetDensityMatrixElement::GetDensityMatrixElement(size_t in_row, size_t in_column):
    m_row(in_row),
    m_column(in_column)
{}

void CalculateDmFidelity::execute(SerializationType* out_result)
{
    assert(out_result != nullptr);
    std::vector<ComplexCoefficient> refVec;
    refVec.reserve(dmElems.size());
    for (const auto& x : dmElems)
    {
        ComplexCoefficient el;
        el.real = x.real();
        el.imag = x.imag();
        refVec.emplace_back(el);
    }

    const double fidelityResult = XACC_QuaC_CalcDensityMatrixFidelity(dmElems.size(), refVec.data());
    SerializationOutputDataType outArchive(*out_result); 
    outArchive(fidelityResult); 
}

CalculateDmFidelity::CalculateDmFidelity(const std::vector<std::complex<double>>& in_refDm):
    dmElems(in_refDm)
{}

void GetDensityMatrixElement::execute(SerializationType* out_result)
{
    assert(out_result != nullptr);
    const auto dmElem = XACC_QuaC_GetDensityMatrixElement(m_row, m_column);
    SerializationOutputDataType outArchive(*out_result); 
    std::complex<double> result { dmElem.real, dmElem.imag };
    outArchive(result); 
}

void FinalizeFunctor::execute(SerializationType* out_result)
{
    XACC_QuaC_Finalize();
}
