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
        const std::vector<double> d_loFreqs { 4.85546, 4.97154 };
        const std::vector<double> u_loFreqs { 4.97154, 4.85546 };
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
            const double ampl = 0.165885;
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
        {
            //  CR drive length = 848*dt
            const int nSamples = 848;
            const double A = 0.47307243770437246;
            const std::string pulseName = "square" + std::to_string(A);
            channelConfigs.addOrReplacePulse(pulseName, QuaC::GaussianSquare(nSamples, A, 32, 720));
            auto cmddef_cnot_01 = provider->createComposite("pulse::cx_0_1");

            auto cr_pulse = std::make_shared<xacc::quantum::Pulse>(pulseName, "u0");
            cr_pulse->setBits({0, 1});
            cr_pulse->setStart(0);
            // Pre- and Post- U3 gates
            // From calibration run:
            // U3(-0.066247, -0.667994, -0.130193), 
            // U3(-1.375548, 1.275995, -1.482965), 
            // U3(0.157918, 0.443962, -0.863998), 
            // U3(0.214932, -0.762514, 0.672792), 

            const auto addDigitalU3 = [&cmddef_cnot_01](size_t qIdx, double in_theta, double in_phi, double in_lambda, int startTime) {
                const auto digitalCmdDefName = "digital::u3_" + std::to_string(qIdx); 

                auto digitalPulse1 = std::make_shared<xacc::quantum::Pulse>(digitalCmdDefName + "_theta");
                auto digitalPulse2 = std::make_shared<xacc::quantum::Pulse>(digitalCmdDefName + "_phi");
                auto digitalPulse3 = std::make_shared<xacc::quantum::Pulse>(digitalCmdDefName + "_lambda");

                digitalPulse1->setBits({qIdx});
                digitalPulse2->setBits({qIdx});
                digitalPulse3->setBits({qIdx});
                digitalPulse1->setStart(startTime);
                digitalPulse2->setStart(startTime);
                digitalPulse3->setStart(startTime);

                xacc::InstructionParameter theta(in_theta);
                digitalPulse1->setParameter(0, theta);

                xacc::InstructionParameter phi(in_phi);
                digitalPulse2->setParameter(0, phi);

                xacc::InstructionParameter lambda(in_lambda);
                digitalPulse3->setParameter(0, lambda);
                cmddef_cnot_01->addInstruction(digitalPulse1);
                cmddef_cnot_01->addInstruction(digitalPulse2);
                cmddef_cnot_01->addInstruction(digitalPulse3);
            };

            // Pre-pulse U3 gates
            // TODO: use actual pulses
            addDigitalU3(0, -0.066247, -0.667994, -0.130193, 0);
            addDigitalU3(1, -1.375548, 1.275995, -1.482965, 0);
            
            // CR pulse
            cmddef_cnot_01->addInstruction(cr_pulse);
            
            // Post-pulse U3 gates
            // TODO: use actual pulses
            addDigitalU3(0, 0.157918, 0.443962, -0.863998, nSamples);
            addDigitalU3(1, 0.214932, -0.762514, 0.672792, nSamples);
            
            xacc::contributeService(cmddef_cnot_01->name(), cmddef_cnot_01);
        }

        setChannelConfigs(channelConfigs);
    }

    const std::string name() const override { return "Default2Q"; }
};
}