#include "PulseSystemModel.hpp"
#include "xacc.hpp"
#include "xacc_service.hpp"
#include "InstructionIterator.hpp"
#include "Pulse.hpp"
#include "json.hpp"

namespace QuaC{

bool HamiltonianModel::fromJson(const std::string& in_jsonString)
{
    auto parser = xacc::getService<QuaC::HamiltonianParsingUtil>("default");
    std::function<void(QuaC::HamiltonianTerm&)> iterFn = [&](QuaC::HamiltonianTerm& in_term) -> void {
       m_terms.emplace_back(in_term.clone());
    };

    // We abort if there is an invalid term
    if (!parser->tryParse(in_jsonString, iterFn))
    {
        m_terms.clear();
        return false;
    }

    auto j = nlohmann::json::parse(in_jsonString);   
    auto qubitDimMap = j["qub"];
    for (auto qbitIter = qubitDimMap.begin(); qbitIter != qubitDimMap.end(); ++qbitIter)
    {
       setQubitDimension(std::stoi(qbitIter.key()), qbitIter.value().get<size_t>());
    }

    return true;
}

size_t HamiltonianModel::getQubitDimension(size_t in_qubitIdx) const 
{
    const auto iter = m_qubitDimension.find(in_qubitIdx);
    if (iter == m_qubitDimension.end() || iter->second < 2)
    {
        return 2;
    }

    return iter->second;
}

bool PulseSystemModel::fromQobjectJson(const std::string& in_jsonString)
{
    auto j = nlohmann::json::parse(in_jsonString);
    const bool isTopLevelConfig = (j.find("hamiltonian") !=  j.end());
    
    auto configJson = isTopLevelConfig ? j : [&](){
        // Note: we will only process one backend per json file
        auto backend = *(j["backends"].begin());
    
        // Get name:
        auto name = backend["name"].get<std::string>();
        if (!name.empty())
        {
            m_name = name;
        }

        // Load Hamiltonain model from Json
        auto config = backend["specificConfiguration"];
        return config;
    }();

    // Load Hamiltonain model from Json
    auto hamiltonian = configJson["hamiltonian"];
    if (!loadHamiltonianJson(hamiltonian.dump()))
    {
        return false;
    }
    
    const double dt =  configJson["dt"].get<double>();
    m_channelConfigs.dt = dt;
    m_channelConfigs.loFregs_dChannels.clear();

    auto loRanges = configJson["qubit_lo_range"];
    for (auto d_iter = loRanges.begin(); d_iter != loRanges.end(); ++d_iter) 
    {
        const double DEFAULT_FREQ = 5.0;
        // TODO: calculate qubit freqs (d channel LO freqs)
        m_channelConfigs.loFregs_dChannels.emplace_back(DEFAULT_FREQ);
    }

    auto uRanges = configJson["u_channel_lo"];
    for (auto u_iter = uRanges.begin(); u_iter != uRanges.end(); ++u_iter) 
    {
        // Adding none freq. but assign the formula appropriately.
        const double DEFAULT_NONE_FREQ = 0.0;
        m_channelConfigs.loFregs_uChannels.emplace_back(DEFAULT_NONE_FREQ);
        uChannelFormula uFreqFormula;
        for (auto it = (*u_iter).begin(); it != (*u_iter).end(); ++it) 
        {
            const int qIdx = (*it)["q"].get<int>();
            const std::vector<double> scale = (*it)["scale"].get<std::vector<double>>();
            assert(scale.size() == 2);
            uFreqFormula.emplace_back(std::make_pair(std::complex<double>(scale[0], scale[1]), qIdx));            
        }
        
        m_channelConfigs.loFregs_uChannelFormulas.emplace_back(std::move(uFreqFormula));
    }

    // Get the pulse library
    auto pulse_library = configJson["defaults"]["pulse_library"];
    
    for (auto pulse_iter = pulse_library.begin(); pulse_iter != pulse_library.end(); ++pulse_iter) 
    {
        auto pulseName = (*pulse_iter)["name"].get<std::string>();
        auto samples = (*pulse_iter)["samples"].get<std::vector<std::vector<double>>>();
        const auto pulseSamples = PulseSamplesToComplexVec(samples);
        if (!pulseSamples.empty())
        {
            m_channelConfigs.addOrReplacePulse(pulseName, std::move(pulseSamples));
        }
    }

    // Import cmd-defs
    if (configJson.find("defaults") !=  configJson.end()) 
    {
        auto cmd_defs = configJson["defaults"]["cmd_def"];
        for (auto cmd_def_iter = cmd_defs.begin(); cmd_def_iter != cmd_defs.end(); ++cmd_def_iter) 
        {
            const auto cmd_def_name = (*cmd_def_iter)["name"].get<std::string>();
            const auto qbits = (*cmd_def_iter)["qubits"].get<std::vector<std::size_t>>();
            std::string tmpName = "pulse::" + cmd_def_name;
            if (cmd_def_name != "measure")
            {
                for (const auto& qb : qbits)
                {
                    tmpName += "_" + std::to_string(qb);
                }
            }

            // Get the composite instruction for the command
            if (xacc::hasContributedService<xacc::Instruction>(tmpName))
            {
                auto cmd_def = xacc::ir::asComposite(xacc::getContributedService<xacc::Instruction>(tmpName));      
                if (!addCommandDef(tmpName, cmd_def))
                {
                    return false;
                }
            }
        }
    }

    return true;
}

bool PulseSystemModel::loadHamiltonianJson(const std::string& in_hamiltonianJsonString)
{
    // Cache the Hamiltonian Json string
    m_hamJson = in_hamiltonianJsonString;
    return m_hamiltonian.fromJson(in_hamiltonianJsonString);
}

void PulseSystemModel::setChannelConfigs(const BackendChannelConfigs& in_config)
{
    m_channelConfigs = in_config;
}

double PulseSystemModel::getQubitT1(size_t in_qubitIdx) const
{
    const auto iter = m_qubitToT1.find(in_qubitIdx);
    if (iter == m_qubitToT1.end())
    {
        return 0.0;
    }

    return iter->second;
}

double PulseSystemModel::getQubitInitialPopulation(size_t in_qubitIdx) const
{
    const auto iter = m_qubitInitialPopulation.find(in_qubitIdx);
    if (iter == m_qubitToT1.end() || iter->second < 0.0 || iter->second > 1.0)
    {
        return 0.0;
    }

    return iter->second;
}

bool PulseSystemModel::addCommandDef(const std::string& in_cmdDefName, const std::shared_ptr<xacc::CompositeInstruction>& in_pulseComposite)
{
    // Check that this composite contains only Pulse instructions
    xacc::InstructionIterator it(in_pulseComposite);
    while (it.hasNext()) 
    {
        auto nextInst = it.next();
        if (nextInst->isEnabled() && !nextInst->isComposite()) 
        {
            auto pulse = std::dynamic_pointer_cast<xacc::quantum::Pulse>(nextInst);
            if (!pulse) 
            {
               return false;
            }

            const auto pulseName = pulse->name();
            if (pulseName != "fc" && pulseName != "acquire" && !m_channelConfigs.hasPulseName(pulseName))
            {
                // Unknown pulse
                return false;
            }
        }
    }
    
    m_pulseCmdDefs.emplace(in_cmdDefName, in_pulseComposite);
    return true;
}

void HamiltonianModel::reset()
{
    m_terms.clear();
    m_qubitDimension.clear();
}

void PulseSystemModel::reset()
{
    m_hamiltonian.reset();
    BackendChannelConfigs defaultConfigs;
    m_channelConfigs = defaultConfigs;
    m_pulseCmdDefs.clear();
    m_qubitToT1.clear();
    m_qubitInitialPopulation.clear();
}
}
