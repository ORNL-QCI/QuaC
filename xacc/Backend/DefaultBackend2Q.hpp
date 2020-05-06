#pragma once

#include "xacc.hpp"
#include "xacc_service.hpp"
#include "PulseSystemModel.hpp"
#include "PulseGen.hpp"
#include "Pulse.hpp"

namespace QuaC {
class Default2Q : public PulseSystemModel
{
public:
    // Default constructor
    Default2Q() 
    {
        // Hamiltonian: full two-qubit transmon Hamiltonian (includes the coupling terms)
        const std::string hamiltonianJson = R"#(
            {
                "description": "Two-qubit Hamiltonian",
                "h_str": ["_SUM[i,0,1,wq{i}*O{i}]", "_SUM[i,0,1,delta{i}*O{i}*(O{i}-I{i})]", "_SUM[i,0,1,omegad{i}*X{i}||D{i}]", "omegad1*X0||U0", "omegad0*X1||U1", "jq0q1*Sp0*Sm1", "jq0q1*Sm0*Sp1"],
                "osc": {},
                "qub": {
                    "0": 2,
                    "1": 2
                },
                "vars": {
                    "wq0": 30.518812656662774, 
                    "wq1": 31.238229295532093,
                    "delta0": -2.011875935,
                    "delta1": -2.008734343,
                    "omegad0": -1.703999855,
                    "omegad1": -1.703999855,
                    "jq0q1": 0.011749557 
                }
            }
        )#";
        loadHamiltonianJson(hamiltonianJson);

        BackendChannelConfigs channelConfigs;
        const std::vector<double> d_loFreqs { 4.857, 4.97154 };
        const std::vector<double> u_loFreqs { 4.972, 4.857 };
        channelConfigs.dt = 0.222;
        channelConfigs.loFregs_dChannels = d_loFreqs;
        channelConfigs.loFregs_uChannels = u_loFreqs;
        auto provider = xacc::getIRProvider("quantum");

        // Pulse library set-up:
        // ======== Single-qubit gate pulse =============     
        // * Q0:
        // X gate:
        {
            auto cmddef_x_0 = provider->createComposite("pulse::x_0");
            cmddef_x_0->setBits({0});
            // These params are pre-calibrated (nl-opt)
            const double ampl = 0.169244;
            const double beta = 1.0;
            const int nSamples = 121;
            const int gamma = 20;

            const std::string pulseName = "drag_" + std::to_string(ampl) + "_" + std::to_string(beta);
            channelConfigs.addOrReplacePulse(pulseName, QuaC::Drag(nSamples, ampl, gamma, beta));
            const std::string channelName = "d0";

            auto pulse = std::make_shared<xacc::quantum::Pulse>(pulseName, channelName);
            pulse->setBits({0});

            cmddef_x_0->addInstructions({ pulse });
            xacc::contributeService(cmddef_x_0->name(), cmddef_x_0);
        }

        // H gate:
        {
            auto cmddef_h_0 = provider->createComposite("pulse::h_0");
            cmddef_h_0->setBits({0});
            // These params are pre-calibrated (nl-opt)
            const double ampl = 0.083113;
            const double beta = 1.0;
            const int nSamples = 121;
            const int gamma = 20;

            const std::string pulseName = "drag_" + std::to_string(ampl) + "_" + std::to_string(beta);
            channelConfigs.addOrReplacePulse(pulseName, QuaC::Drag(nSamples, ampl, gamma, beta));
            const std::string channelName = "d0";

            auto pulse1 = std::make_shared<xacc::quantum::Pulse>(pulseName, channelName);
            pulse1->setBits({0});

            xacc::InstructionParameter fcParameter(M_PI_2);
            auto fcInst1 = std::make_shared<xacc::quantum::Pulse>("fc", channelName);
            fcInst1->setParameter(0, fcParameter);
            fcInst1->setBits({0});
            fcInst1->setStart(nSamples);

            auto fcInst2 = std::make_shared<xacc::quantum::Pulse>("fc", "u0");
            fcInst2->setParameter(0, fcParameter);
            fcInst2->setBits({0});
            fcInst2->setStart(nSamples);
            
            auto pulse2 = std::make_shared<xacc::quantum::Pulse>(pulseName, channelName);
            pulse2->setBits({0});
            pulse2->setStart(nSamples);

            cmddef_h_0->addInstructions({ pulse1, fcInst1, fcInst2, pulse2 });
            xacc::contributeService(cmddef_h_0->name(), cmddef_h_0);
        }

        // U3 gate
        {
            // These params are pre-calibrated (nl-opt)
            const double ampl = 0.083113;
            const double beta = 1.0;
            const int nSamples = 121;
            const int gamma = 20;

            const std::string pulseName = "drag_" + std::to_string(ampl) + "_" + std::to_string(beta);
            channelConfigs.addOrReplacePulse(pulseName, QuaC::Drag(nSamples, ampl, gamma, beta));
            
            auto cmddef_u3_0 = provider->createComposite("pulse::u3_0");
            cmddef_u3_0->addVariables({"P0", "P1", "P2"});
            cmddef_u3_0->setBits({0});

            auto fc_P0_d = std::make_shared<xacc::quantum::Pulse>("fc", "d0");
            xacc::InstructionParameter fcParameterP0("pi - P0");
            fc_P0_d->setParameter(0, fcParameterP0);
            fc_P0_d->setStart(nSamples);
            fc_P0_d->setBits({0});

            auto fc_P0_u = std::make_shared<xacc::quantum::Pulse>("fc", "u0");
            fc_P0_u->setParameter(0, fcParameterP0);
            fc_P0_u->setStart(nSamples);
            fc_P0_u->setBits({0});


            auto fc_P1_d = std::make_shared<xacc::quantum::Pulse>("fc", "d0");
            xacc::InstructionParameter fcParameterP1("P1 - pi/2");
            fc_P1_d->setParameter(0, fcParameterP1);
            fc_P1_d->setStart(nSamples + nSamples);
            fc_P1_d->setBits({0});

            auto fc_P1_u = std::make_shared<xacc::quantum::Pulse>("fc", "u0");
            fc_P1_u->setParameter(0, fcParameterP1);
            fc_P1_u->setStart(nSamples + nSamples);
            fc_P1_u->setBits({0});
            
            auto fc_P2_d = std::make_shared<xacc::quantum::Pulse>("fc", "d0");
            xacc::InstructionParameter fcParameterP2("P2 - pi/2");
            fc_P2_d->setParameter(0, fcParameterP2);
            fc_P2_d->setStart(0);
            fc_P2_d->setBits({0});

            auto fc_P2_u = std::make_shared<xacc::quantum::Pulse>("fc", "u0");
            fc_P2_u->setParameter(0, fcParameterP2);
            fc_P2_u->setStart(0);
            fc_P2_u->setBits({0});

            auto firstPulseX_Pi_2 = std::make_shared<xacc::quantum::Pulse>(pulseName, "d0");
            firstPulseX_Pi_2->setStart(0);
            firstPulseX_Pi_2->setBits({0});
            
            auto secondPulseX_Pi_2 = std::make_shared<xacc::quantum::Pulse>(pulseName, "d0");
            secondPulseX_Pi_2->setStart(nSamples);
            secondPulseX_Pi_2->setBits({0});

            // Add pulse sequence to the cmd-def
            cmddef_u3_0->addInstructions({fc_P2_d, fc_P2_u, firstPulseX_Pi_2, fc_P0_d, fc_P0_u, secondPulseX_Pi_2, fc_P1_d, fc_P1_u });
            xacc::contributeService(cmddef_u3_0->name(), cmddef_u3_0);
        }
        // * Q1:
        // X gate:
        {
            auto cmddef_x_1 = provider->createComposite("pulse::x_1");
            cmddef_x_1->setBits({1});
            // These params are pre-calibrated (nl-opt)
            const double ampl = 0.169244;
            const double beta = 1.0;
            const int nSamples = 121;
            const int gamma = 20;

            const std::string pulseName = "drag_" + std::to_string(ampl) + "_" + std::to_string(beta);
            channelConfigs.addOrReplacePulse(pulseName, QuaC::Drag(nSamples, ampl, gamma, beta));
            const std::string channelName = "d1";

            auto pulse = std::make_shared<xacc::quantum::Pulse>(pulseName, channelName);
            pulse->setBits({1});

            cmddef_x_1->addInstructions({ pulse });
            xacc::contributeService(cmddef_x_1->name(), cmddef_x_1);
        }

        // H gate:
        {
            auto cmddef_h_1 = provider->createComposite("pulse::h_1");
            cmddef_h_1->setBits({1});
            // These params are pre-calibrated (nl-opt)
            const double ampl = 0.0813729;
            const double beta = 1.0;
            const int nSamples = 121;
            const int gamma = 20;

            const std::string pulseName = "drag_" + std::to_string(ampl) + "_" + std::to_string(beta);
            channelConfigs.addOrReplacePulse(pulseName, QuaC::Drag(nSamples, ampl, gamma, beta));
            const std::string channelName = "d1";

            auto pulse1 = std::make_shared<xacc::quantum::Pulse>(pulseName, channelName);
            pulse1->setBits({1});

            xacc::InstructionParameter fcParameter(M_PI_2);
            auto fcInst1 = std::make_shared<xacc::quantum::Pulse>("fc", channelName);
            fcInst1->setParameter(0, fcParameter);
            fcInst1->setBits({1});
            fcInst1->setStart(nSamples);

            auto fcInst2 = std::make_shared<xacc::quantum::Pulse>("fc", "u1");
            fcInst2->setParameter(0, fcParameter);
            fcInst2->setBits({1});
            fcInst2->setStart(nSamples);
            
            auto pulse2 = std::make_shared<xacc::quantum::Pulse>(pulseName, channelName);
            pulse2->setBits({1});
            pulse2->setStart(nSamples);

            cmddef_h_1->addInstructions({ pulse1, fcInst1, fcInst2, pulse2 });
            xacc::contributeService(cmddef_h_1->name(), cmddef_h_1);
        }

        // U3 gate
        {
            // These params are pre-calibrated (nl-opt)
            const double ampl = 0.0813729;
            const double beta = 1.0;
            const int nSamples = 121;
            const int gamma = 20;

            const std::string pulseName = "drag_" + std::to_string(ampl) + "_" + std::to_string(beta);
            channelConfigs.addOrReplacePulse(pulseName, QuaC::Drag(nSamples, ampl, gamma, beta));
            
            auto cmddef_u3_1 = provider->createComposite("pulse::u3_1");
            cmddef_u3_1->addVariables({"P0", "P1", "P2"});
            cmddef_u3_1->setBits({1});

            auto fc_P0_d = std::make_shared<xacc::quantum::Pulse>("fc", "d1");
            xacc::InstructionParameter fcParameterP0("pi - P0");
            fc_P0_d->setParameter(0, fcParameterP0);
            fc_P0_d->setStart(nSamples);
            fc_P0_d->setBits({1});

            auto fc_P0_u = std::make_shared<xacc::quantum::Pulse>("fc", "u1");
            fc_P0_u->setParameter(0, fcParameterP0);
            fc_P0_u->setStart(nSamples);
            fc_P0_u->setBits({1});


            auto fc_P1_d = std::make_shared<xacc::quantum::Pulse>("fc", "d1");
            xacc::InstructionParameter fcParameterP1("P1 - pi/2");
            fc_P1_d->setParameter(0, fcParameterP1);
            fc_P1_d->setStart(nSamples + nSamples);
            fc_P1_d->setBits({1});

            auto fc_P1_u = std::make_shared<xacc::quantum::Pulse>("fc", "u1");
            fc_P1_u->setParameter(0, fcParameterP1);
            fc_P1_u->setStart(nSamples + nSamples);
            fc_P1_u->setBits({1});
            
            auto fc_P2_d = std::make_shared<xacc::quantum::Pulse>("fc", "d1");
            xacc::InstructionParameter fcParameterP2("P2 - pi/2");
            fc_P2_d->setParameter(0, fcParameterP2);
            fc_P2_d->setStart(0);
            fc_P2_d->setBits({1});

            auto fc_P2_u = std::make_shared<xacc::quantum::Pulse>("fc", "u1");
            fc_P2_u->setParameter(0, fcParameterP2);
            fc_P2_u->setStart(0);
            fc_P2_u->setBits({1});

            auto firstPulseX_Pi_2 = std::make_shared<xacc::quantum::Pulse>(pulseName, "d1");
            firstPulseX_Pi_2->setStart(0);
            firstPulseX_Pi_2->setBits({1});
            
            auto secondPulseX_Pi_2 = std::make_shared<xacc::quantum::Pulse>(pulseName, "d1");
            secondPulseX_Pi_2->setStart(nSamples);
            secondPulseX_Pi_2->setBits({1});

            // Add pulse sequence to the cmd-def
            cmddef_u3_1->addInstructions({fc_P2_d, fc_P2_u, firstPulseX_Pi_2, fc_P0_d, fc_P0_u, secondPulseX_Pi_2, fc_P1_d, fc_P1_u });
            xacc::contributeService(cmddef_u3_1->name(), cmddef_u3_1);
        }


        // CNOT gate: 
        // NOTE: We temporary use digital U3 gates/pseudo pulses to
        // correct local roration errors pre- and post- CR pulses.
        // We will replace them  which actual pulses later.
        // Rationale: This will save significant simulation time during development/testing.
        // ===============================================     


        setChannelConfigs(channelConfigs);
    }

    const std::string name() const override { return "Default2Q"; }
};
}