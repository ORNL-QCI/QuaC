#include "QuaC_Accelerator.hpp"
#include "QuaC_Pulse_Visitor.hpp"
#include "Pulse.hpp"
#include "PulseSystemModel.hpp"
#include "BackendRegistry.hpp"
#include "xacc_service.hpp"

namespace QuaC {
    void QuaC_Accelerator::initialize(const HeterogeneousMap& params)  
    {        
        if (params.stringExists("backend"))
        {
            const auto backendName = params.getString("backend");
            std::cout << "Backend name: " << backendName << "\n";
            // If the name starts with "ibmq", e.g., "ibmq_armonk",
            // delagate to the IBM remote accelerator to initialize pulse library.
            if (backendName.size() > 4 && backendName.substr(0, 4) == "ibmq") {
                auto ibmAcc = xacc::getAccelerator("ibm:" + backendName);
                ibmAcc->contributeInstructions();

                m_systemModel = std::make_shared<PulseSystemModel>();
                const auto backendProps = ibmAcc->getProperties().get<std::string>("total-json");
                std::cout << "Backend JSON:\n" << backendProps << "\n";
                if (!m_systemModel->fromQobjectJson(backendProps))
                {
                    xacc::error("Failed to initialize pulse system model from the JSON file.");
                    return;
                }
            }
            else 
            {
                // Not an IBM backend: try to look up from our sample backend models.
                auto backendRegistry = xacc::getService<BackendRegistry>("default");
                m_systemModel = backendRegistry->getSystemModel(backendName);
                if (!m_systemModel)
                {
                    xacc::error("Invalid backend named '" + backendName + "' was requested.");
                }

                if (params.keyExists<std::vector<double>>("initial-population"))
                {
                    const std::vector<double> initialPops = params.get<std::vector<double>>("initial-population");
                    for (size_t i = 0; i < initialPops.size(); ++i)
                    {
                        m_systemModel->setQubitInitialPopulation(i, initialPops[i]);
                    }
                }
            }
        }
        else
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
        }
        
        assert(m_systemModel);
        m_params = params;       
        m_pulseVisitor = std::make_shared<PulseVisitor>();
    }

    void QuaC_Accelerator::execute(std::shared_ptr<AcceleratorBuffer> buffer, const std::shared_ptr<CompositeInstruction> compositeInstruction)  
    {
        // Handle pulse-level IR transformation:
        // In this mode: we will return the Hamiltonian in the buffer rather than executing the circuit.
        // TODO: formalize this key
        if (buffer->hasExtraInfoKey("ir-transform") && compositeInstruction == nullptr)
        {
            m_pulseVisitor->retrievePulseSystemModel(buffer, m_systemModel.get(), m_params);
            return;
        }
        
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
        
        // If there is a target density matrix, 
        // calculate the fidelity b/w the result Dm and target dm.
        if (buffer->hasExtraInfoKey("target-dm-real") && buffer->hasExtraInfoKey("target-dm-imag"))
        {
            const std::vector<double> targetDmReal = (*buffer)["target-dm-real"].as<std::vector<double>>();
            const std::vector<double> targetDmImag = (*buffer)["target-dm-imag"].as<std::vector<double>>();

            assert(targetDmReal.size() == targetDmImag.size());
            std::vector<std::complex<double>> targetDm;
            targetDm.reserve(targetDmReal.size());
            for (size_t i = 0; i < targetDmReal.size(); ++i)
            {
                targetDm.emplace_back(std::complex<double> { targetDmReal[i], targetDmImag[i] });
            }

            const auto fidelity = m_pulseVisitor->calcFidelity(targetDm);
            buffer->addExtraInfo("fidelity", fidelity);
        }

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
        
        const auto keyExists = [](const nlohmann::json& in_json, const std::string& in_key) {
            return in_json.find(in_key) != in_json.end();
        };

        // Note: There are currently two types of Json files that can be used to contribute pulse instructions here
        // (1) A full pulse backend config file which includes the system dynamical model (aka Hamiltonian) and the *pulse library*
        // (2) A subset of (1) which only contains the pulse library.
        // In this contributeInstructions() method, we only need (2) but just in case the input Json file is of type (1),
        // just look up the pulse library field accordingly.
        
        // Always add frame change and acquire instructions
        auto fc = std::make_shared<Pulse>("fc");
        xacc::contributeService("fc", fc);
        auto aq = std::make_shared<Pulse>("acquire");
        xacc::contributeService("acquire", aq);

        // This is a full config Json file:
        if (keyExists(j, "backends"))
        {
            auto backends = j["backends"];      
            for (auto it = backends.begin(); it != backends.end(); ++it) 
            {
                // Get the pulse library
                auto pulse_library = (*it)["specificConfiguration"]["defaults"]["pulse_library"];
                contributePulseInstructions(pulse_library);

                // Import command defs (sequence of pulses)
                auto cmd_defs = (*it)["specificConfiguration"]["defaults"]["cmd_def"];
                contributeCmdDefInstructions(cmd_defs);
            }
        }
        else
        {
            // The Json is a pulse library only file
            if (keyExists(j, "pulse_library"))
            {
                contributePulseInstructions(j["pulse_library"]);
            }
            if (keyExists(j, "cmd_def"))
            {
                contributeCmdDefInstructions(j["cmd_def"]);
            }
        }    
    }

    void QuaC_Accelerator::contributePulseInstructions(const nlohmann::json& in_pulseLibJson)
    {
        if (!in_pulseLibJson.is_array())
        {
            xacc::error("Pulse library must be an array object.\n");
            return;
        }

        for (auto pulse_iter = in_pulseLibJson.begin(); pulse_iter != in_pulseLibJson.end(); ++pulse_iter) 
        {
            auto pulse_name = (*pulse_iter)["name"].get<std::string>();
            auto samples = (*pulse_iter)["samples"].get<std::vector<std::vector<double>>>();
            auto pulse = std::make_shared<xacc::quantum::Pulse>(pulse_name);
            pulse->setSamples(samples);
            xacc::contributeService(pulse_name, pulse);
        }
    }

    void QuaC_Accelerator::contributeCmdDefInstructions(const nlohmann::json& in_cmdDefJson)
    {
        if (!in_cmdDefJson.is_array())
        {
            xacc::error("Pulse cmd-def list must be an array object.\n");
            return;
        }

        auto provider = xacc::getIRProvider("quantum");
        for (auto cmd_def_iter = in_cmdDefJson.begin(); cmd_def_iter != in_cmdDefJson.end(); ++cmd_def_iter) 
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

