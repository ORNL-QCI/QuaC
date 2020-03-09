#pragma once
#include "Serialization.hpp"
#include "PulseChannelController.hpp"

struct FunctorBase
{
    virtual void execute(SerializationType* out_result = nullptr) = 0;
    virtual std::string name() const = 0;
    virtual ~FunctorBase() {}
};

struct InitializeFunctor: public FunctorBase
{
    virtual void execute(SerializationType* out_result = nullptr) override;
    virtual std::string name() const override { return "InitializeFunctor"; }

    InitializeFunctor() {}
    InitializeFunctor(int in_nbQubit, 
                        const std::vector<int>& in_qbitDims, 
                        const std::unordered_map<int, double>& in_qbitDecays = {}, 
                        const std::unordered_map<int, double>& in_qbitInitialPopulations = {},
                        bool in_verbose = false,
                        double in_monitorDt = -1.0);
    
    template<class Archive>
    void serialize(Archive& archive)
    {
        archive(m_nbQubit, m_qbitDims, m_qbitDecays, m_qbitInitialPops, m_verbose, m_monitorDt); 
    }
    int m_nbQubit;
    std::vector<int> m_qbitDims;
    std::vector<double> m_qbitDecays;
    std::vector<double> m_qbitInitialPops;
    bool m_verbose;
    double m_monitorDt;
};
DECLARE_CORE_TYPE(FunctorBase, InitializeFunctor)

struct AddHamiltonianTerm: public FunctorBase
{
    virtual void execute(SerializationType* out_result = nullptr) override;
    virtual std::string name() const override { return "AddHamiltonianTerm"; }

    AddHamiltonianTerm() {}
    AddHamiltonianTerm(const std::complex<double>& in_coeff, const std::vector<std::pair<std::string, int>>& in_ops, int in_channelId = -1);
    
    template<class Archive>
    void serialize(Archive& archive)
    {
        archive(m_channelId, m_coeff, m_ops); 
    }
    int m_channelId;
    std::complex<double> m_coeff;
    std::vector<std::pair<std::string, int>> m_ops;
}; 
DECLARE_CORE_TYPE(FunctorBase, AddHamiltonianTerm)


struct AddGateU3: public FunctorBase
{
    virtual void execute(SerializationType* out_result = nullptr) override;
    virtual std::string name() const override { return "AddGateU3"; }
    
    AddGateU3() {}
    AddGateU3(int in_qubitIdx, double in_theta, double in_phi, double in_lambda, double in_startTime);
    
    template<class Archive>
    void serialize(Archive& archive)
    {
        archive(m_qubitIdx, m_theta, m_phi, m_lambda, m_startTime); 
    }
    int m_qubitIdx;
    double m_theta; 
    double m_phi;
    double m_lambda;
    double m_startTime;
}; 
DECLARE_CORE_TYPE(FunctorBase, AddGateU3)

struct AddGateCNOT: public FunctorBase
{
    virtual void execute(SerializationType* out_result = nullptr) override;
    virtual std::string name() const override { return "AddGateCNOT"; }
    
    AddGateCNOT() {}
    AddGateCNOT(int in_ctrlIdx, int in_targetIdx, double in_startTime);
    
    template<class Archive>
    void serialize(Archive& archive)
    {
        archive(m_ctrlIdx, m_targetIdx, m_startTime); 
    }

    int m_ctrlIdx;
    int m_targetIdx;
    double m_startTime;
}; 
DECLARE_CORE_TYPE(FunctorBase, AddGateCNOT)

struct StartTimestepping: public FunctorBase
{
    virtual void execute(SerializationType* out_result = nullptr) override;
    virtual std::string name() const override { return "StartTimestepping"; }

    StartTimestepping(const PulseChannelController& in_pulseDataProvider, double in_dt, double in_stopTime, int in_stepMax, bool in_adaptive = true);
    StartTimestepping() {}
    
    template<class Archive>
    void serialize(Archive& archive)
    {
        archive(m_pulseDataProvider, m_dt, m_stopTime, m_stepMax, m_adaptive); 
    }
    PulseChannelController m_pulseDataProvider;
    double m_dt;
    double m_stopTime;
    int m_stepMax;
    bool m_adaptive;
}; 
DECLARE_CORE_TYPE(FunctorBase, StartTimestepping)

struct CalculateBipartiteConcurrence: public FunctorBase
{
    size_t m_qubit1;
    size_t m_qubit2;

    virtual void execute(SerializationType* out_result = nullptr) override;
    virtual std::string name() const override { return "CalculateBipartiteConcurrence"; }

    CalculateBipartiteConcurrence(){};
    CalculateBipartiteConcurrence(size_t in_q1, size_t in_q2);

    template<class Archive>
    void serialize(Archive& archive) 
    {
        archive(m_qubit1, m_qubit2); 
    }
};
DECLARE_CORE_TYPE(FunctorBase, CalculateBipartiteConcurrence)


struct GetDensityMatrixElement: public FunctorBase
{
    size_t m_row;
    size_t m_column;

    virtual std::string name() const override { return "GetDensityMatrixElement"; }
    virtual void execute(SerializationType* out_result = nullptr) override;
    GetDensityMatrixElement(){};
    GetDensityMatrixElement(size_t in_row, size_t in_column);

    template<class Archive>
    void serialize(Archive& archive) 
    {
        archive(m_row, m_column); 
    }
};
DECLARE_CORE_TYPE(FunctorBase, GetDensityMatrixElement)

struct FinalizeFunctor: public FunctorBase
{
    virtual std::string name() const override { return "FinalizeFunctor"; }
    virtual void execute(SerializationType* out_result = nullptr) override;
    FinalizeFunctor(){};
    template<class Archive>
    void serialize(Archive& archive) {}
};
DECLARE_CORE_TYPE(FunctorBase, FinalizeFunctor)

struct TimeSteppingData {
    double time;
    std::vector<double> channelData;
    std::vector<double> populations;
    std::vector<double> pauliExpectations;

    template<class Archive>
    void serialize(Archive& archive)
    {
        archive(time, channelData, populations, pauliExpectations); 
    }
};

struct SimResult
{
    std::vector<double> finalPopulations;
    std::vector<TimeSteppingData> tsData; 
    std::vector<std::complex<double>> dmDiagElems;
    // If not empty, absolute path to the csv file
    // that contains time-stepping data
    std::string csvFileName;
    template<class Archive>
    void serialize(Archive& archive)
    {
        archive(finalPopulations, tsData, dmDiagElems, csvFileName); 
    }
};
