#include "QuaC_Accelerator.hpp"
#include "QuaC_Pulse_Visitor.hpp"
#include "json.hpp"
#include "Pulse.hpp"
#include "PulseSystemModel.hpp"

namespace QuaC {
    void QuaC_Accelerator::initialize(const HeterogeneousMap& params)  
    {        
        if (params.stringExists("config-json-path"))
        {
            const auto jsonFileName = params.getString("config-json-path");
            std::ifstream backendFile(jsonFileName);
            std::string jjson((std::istreambuf_iterator<char>(backendFile)), std::istreambuf_iterator<char>());
            // Contribute all the pulse instruction defined in the config file.
            contributeInstructions(jjson);
            m_systemModel = std::make_shared<PulseSystemModel>();
            if (!m_systemModel->fromQobjectJson(jjson))
            {
                xacc::error("Failed to initialize pulse system model from the JSON file.");
                return;
            }
        }
        else if (params.pointerLikeExists<PulseSystemModel>("system-model")) 
        {
            PulseSystemModel* systemModel = params.getPointerLike<PulseSystemModel>("system-model");
            // we don't own this one, don't try to delete it.
            m_systemModel = std::shared_ptr<PulseSystemModel>(systemModel , [](PulseSystemModel* ptr){});
        } 
        else if (params.stringExists("system-model"))
        {
            // Request the system model as a service
            // std::cout << "Using system model: " << params.getString("system-model") << "\n";
            auto serviceRef = xacc::getService<PulseSystemModel>(params.getString("system-model"));
            m_systemModel = std::shared_ptr<PulseSystemModel>(serviceRef.get(), [](PulseSystemModel* ptr){}); 
        }
        else
        {
            xacc::error("Either a PulseSystemModel object or a path to the JSON configs file is required.");
            return;
        }
        
        assert(m_systemModel);
        m_params = params;       
        m_pulseVisitor = std::make_shared<PulseVisitor>();
    }

    void QuaC_Accelerator::execute(std::shared_ptr<AcceleratorBuffer> buffer, const std::shared_ptr<CompositeInstruction> compositeInstruction)  
    {
        m_pulseVisitor->initialize(buffer, m_systemModel.get(), m_params);
        // Walk the IR tree, and visit each node
        InstructionIterator it(compositeInstruction);
        while (it.hasNext()) 
        {
            auto nextInst = it.next();
            if (nextInst->isEnabled()) 
            {
                nextInst->accept(m_pulseVisitor);
            }
        }            
        m_pulseVisitor->solve();
        m_pulseVisitor->finalize();
    }

    void QuaC_Accelerator::execute(std::shared_ptr<AcceleratorBuffer> buffer, const std::vector<std::shared_ptr<CompositeInstruction>> compositeInstructions)  
    {  
       for (const auto& inst: compositeInstructions)
       {
            auto tmpBuffer = std::make_shared<xacc::AcceleratorBuffer>(inst->name(), buffer->size());
            execute(tmpBuffer, inst);
            buffer->appendChild(inst->name(), tmpBuffer);
        }
    }

    void QuaC_Accelerator::contributeInstructions(const std::string& custom_json_config) 
    {
        if (custom_json_config.empty())
        {
            return;
        }

        auto provider = xacc::getIRProvider("quantum");
        auto j = nlohmann::json::parse(custom_json_config);
        auto backends = j["backends"];

        // Add frame change and acquire
        auto fc = std::make_shared<Pulse>("fc");
        xacc::contributeService("fc", fc);
        auto aq = std::make_shared<Pulse>("acquire");
        xacc::contributeService("acquire", aq);

        for (auto it = backends.begin(); it != backends.end(); ++it) 
        {
            // Get the pulse library
            auto pulse_library = (*it)["specificConfiguration"]["defaults"]["pulse_library"];
            int counter = 0;
            for (auto pulse_iter = pulse_library.begin(); pulse_iter != pulse_library.end(); ++pulse_iter) 
            {
                auto pulse_name = (*pulse_iter)["name"].get<std::string>();
                auto samples = (*pulse_iter)["samples"].get<std::vector<std::vector<double>>>();

                // std::cout << counter << ", Pulse: " << pulse_name << " number of samples = " << samples.size() << " \n";

                auto pulse = std::make_shared<xacc::quantum::Pulse>(pulse_name);
                pulse->setSamples(samples);
                xacc::contributeService(pulse_name, pulse);
                counter++;
            }

            // Import command defs (sequence of pulses)
            auto cmd_defs = (*it)["specificConfiguration"]["defaults"]["cmd_def"];
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
                // Create a composite instruction for the command
                auto cmd_def = provider->createComposite(tmpName);
                // Add params if they are parameterized commands
                if (cmd_def_name == "u3") 
                {
                    cmd_def->addVariables({"P0", "P1", "P2"});
                } 
                else if (cmd_def_name == "u1") 
                {
                    cmd_def->addVariables({"P0"});
                } 
                else if (cmd_def_name == "u2") 
                {
                    cmd_def->addVariables({"P0", "P1"});
                }

                auto sequence = (*cmd_def_iter)["sequence"];
                for (auto seq_iter = sequence.begin(); seq_iter != sequence.end(); ++seq_iter) 
                {
                    const auto inst_name = (*seq_iter)["name"].get<std::string>();
                    auto inst = xacc::getContributedService<Instruction>(inst_name);

                    if (inst_name != "acquire") 
                    {
                        const auto channel = (*seq_iter)["ch"].get<std::string>();
                        const auto t0 = (*seq_iter)["t0"].get<int>();
                        inst->setBits(qbits);
                        inst->setChannel(channel);
                        inst->setStart(t0);

                        if ((*seq_iter).find("phase") != (*seq_iter).end()) 
                        {
                            // we have phase too
                            auto p = (*seq_iter)["phase"];
                            if (p.is_string()) 
                            {
                                // this is a variable we have to keep track of
                                auto ptmp = p.get<std::string>();
                                // get true variable
                                ptmp.erase(std::remove_if(ptmp.begin(), ptmp.end(), [](char ch) { return ch == '(' || ch == ')'; }), ptmp.end());
                                InstructionParameter phase(ptmp);
                                inst->setParameter(0, phase);
                            } 
                            else 
                            {
                                InstructionParameter phase(p.get<double>());
                                inst->setParameter(0, phase);
                            }
                        }
                    }
                    cmd_def->addInstruction(inst);
                }
                cmd_def->setBits(qbits);

                xacc::contributeService(tmpName, cmd_def);
            }
        }
    }
}

