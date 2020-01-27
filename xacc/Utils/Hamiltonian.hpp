#pragma once

#include <string>
#include <complex>
#include <vector>
#include <memory>
#include <functional>
#include <unordered_map>
#include "Identifiable.hpp"

namespace QuaC {
enum class Operator { X, Y, Z, SP, SM, I, O, P, NA };

inline std::string OperatorToString(Operator op)
{
    switch(op)
    {
        case Operator::X: return "X";
        case Operator::Y: return "Y";
        case Operator::Z: return "Z";
        case Operator::SP: return "SP";
        case Operator::SM: return "SM";
        case Operator::I: return "I";
        case Operator::O: return "O";
        case Operator::P: return "P";
        default: return "";
    }
}

inline Operator ConvertOperatorFromString(const std::string& str)
{
    static std::unordered_map<std::string, Operator> strToOpMap;
    if (strToOpMap.empty())
    {
        for (int enumIter = static_cast<int>(Operator::X); enumIter != static_cast<int>(Operator::NA); ++enumIter)
        {
            auto opEnum = static_cast<Operator>(enumIter);
            strToOpMap.emplace(OperatorToString(opEnum), opEnum);
        }
    }
    
    const auto iter = strToOpMap.find(str);

    if (iter != strToOpMap.end())
    {
        return iter->second;
    }
    else
    {
        return Operator::NA;
    }    
}

// A qubit operator is a pair of operator type and qubit index
using QubitOp = std::pair<Operator, size_t>;
// Variable maps (i.e. to resolve variables in the Hamiltonian) 
using VarsMap = std::unordered_map<std::string, double>;

struct IChannelNameResolver
{
    virtual int GetChannelId(const std::string& in_channelName) = 0;
};

class HamiltonianTerm 
{
public:
    // Apply the Hamiltonian term to the backend.
    virtual void apply(IChannelNameResolver* in_channelResolver) = 0;
    virtual ~HamiltonianTerm() {}
};

class HamiltonianSumTerm: public HamiltonianTerm
{
public:
    static std::unique_ptr<HamiltonianSumTerm> fromString(const std::string& in_string, const VarsMap& in_vars);

    HamiltonianSumTerm(std::vector<std::unique_ptr<HamiltonianTerm>>&& terms, bool isTimeDependent):
        m_terms(std::move(terms)), 
        m_isTimeDependent(isTimeDependent)
    {}

    virtual void apply(IChannelNameResolver* in_channelResolver) override;

private:
    std::vector<std::unique_ptr<HamiltonianTerm>> m_terms;
    bool m_isTimeDependent;
};

class HamiltonianTimeIndependentTerm: public HamiltonianTerm
{
public:
    static std::unique_ptr<HamiltonianTimeIndependentTerm> fromString(const std::string& in_string, const VarsMap& in_vars);

    HamiltonianTimeIndependentTerm(const std::complex<double>& in_coeff, const std::vector<QubitOp>& in_ops):
        m_coefficient(in_coeff),
        m_operators(in_ops)
    {}
    
    virtual void apply(IChannelNameResolver* in_channelResolver) override;

private:
    std::complex<double> m_coefficient;
    std::vector<QubitOp> m_operators;
};

class HamiltonianTimeDependentTerm: public HamiltonianTerm
{
public:    
    // Format <op>||<ch> (channel is Di or Ui)
    static std::unique_ptr<HamiltonianTimeDependentTerm> fromString(const std::string& in_string, const VarsMap& in_vars);
    
    HamiltonianTimeDependentTerm(const std::string& in_channelName, double in_coeff, const QubitOp& in_op):
        m_channelName(in_channelName),
        m_coefficient(in_coeff),
        m_operator(in_op)
    {}

    virtual void apply(IChannelNameResolver* in_channelResolver) override;

private:
    std::string m_channelName;
    double m_coefficient;
    // Note: currently, we only support single operator in time-dependent terms
    QubitOp m_operator;
};

// Hamiltonian Parsing utility for IBM Open Pulse format
class HamiltonianParsingUtil : public xacc::Identifiable 
{
public:
  // Null if cannot parse
  std::unique_ptr<HamiltonianTerm> tryParse(const std::string& in_expr, const VarsMap& in_vars);
  bool tryParse(const std::string& in_jsonString, std::function<void(HamiltonianTerm&)> in_forEachTermFn);
  const std::string name() const override { return "default"; }
  const std::string description() const override { return "Parser for Open Pulse Hamiltonian terms"; }
};
}