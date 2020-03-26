#pragma once

#include "Hamiltonian.hpp"
#include "PulseChannelController.hpp"
#include "CompositeInstruction.hpp"

namespace QuaC{
// A Hamiltonian model contains Hamiltonian data for pulse simulation
class HamiltonianModel
{
public:
    bool fromJson(const std::string& in_jsonString);

    std::vector<std::unique_ptr<HamiltonianTerm>>& getTerms() { return m_terms; }
    void setQubitDimension(size_t in_qubitIdx, size_t in_dimension) { m_qubitDimension[in_qubitIdx] = in_dimension; }
    size_t getQubitDimension(size_t in_qubitIdx) const;
    void reset();
private:
    std::vector<std::unique_ptr<HamiltonianTerm>> m_terms;
    std::unordered_map<size_t, size_t> m_qubitDimension;
};


// A PulseSystemModel encapsulates all model parameters necessary for dynamical simulation.
class PulseSystemModel : public xacc::Identifiable 
{
public:
    PulseSystemModel(const std::string& in_name = "default") :
        m_name(in_name)
    {}

    const std::string name() const override { return m_name; }
    const std::string description() const override { return "Pulse-level System Model"; }
    
    // Try to load the entire system model from Json (must be in the IBM Open Pulse format)
    bool fromQobjectJson(const std::string& in_jsonString);
    
    // Try to load the Hamiltonian Json to initialize the Hamiltonian model.
    bool loadHamiltonianJson(const std::string& in_hamiltonianJsonString);

    // Retrieve the OpenPulse Hamiltonian Json
    std::string getHamiltonianJson() const { return m_hamJson; }
    
    // Set the channel configs 
    void setChannelConfigs(const BackendChannelConfigs& in_config);
    
    // Add a command def (as a composite instruction consisting of pulse instructions)
    // Return false if failed, e.g. the pulse instruction in the composite is not a valid one.
    bool addCommandDef(const std::string& in_cmdDefName, const std::shared_ptr<xacc::CompositeInstruction>& in_pulseComposite);

    BackendChannelConfigs& getChannelConfigs() { return m_channelConfigs; } 
    
    HamiltonianModel& getHamiltonian() { return m_hamiltonian; }
    
    // Set T1 on a qubit (this info can also be parsed from the Json)
    void setQubitT1(size_t in_qubitIdx, double in_qubitT1) { m_qubitToT1[in_qubitIdx] = in_qubitT1; }
    double getQubitT1(size_t in_qubitIdx) const;

    void setQubitInitialPopulation(size_t in_qubitIdx, double in_initialPopulation) { m_qubitInitialPopulation[in_qubitIdx] = in_initialPopulation; }
    double getQubitInitialPopulation(size_t in_qubitIdx) const;
    void reset();
private:
    std::string m_name;
    HamiltonianModel m_hamiltonian;
    BackendChannelConfigs m_channelConfigs;
    std::unordered_map<std::string, std::shared_ptr<xacc::CompositeInstruction>> m_pulseCmdDefs;
    std::unordered_map<size_t, double> m_qubitToT1;
    std::unordered_map<size_t, double> m_qubitInitialPopulation;
    std::string m_hamJson;
};
}
