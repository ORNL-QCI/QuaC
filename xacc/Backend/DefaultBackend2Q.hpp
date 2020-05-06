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
        const std::vector<double> d_loFreqs { 4.857, 4.972 };
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
            const double fcAngle = 1.56906;
            const int nSamples = 121;
            const int gamma = 20;

            const std::string pulseName = "drag_" + std::to_string(ampl) + "_" + std::to_string(beta);
            channelConfigs.addOrReplacePulse(pulseName, QuaC::Drag(nSamples, ampl, gamma, beta));
            const std::string channelName = "d0";

            auto pulse = std::make_shared<xacc::quantum::Pulse>(pulseName, channelName);
            pulse->setBits({0});

            auto fcPre = std::make_shared<xacc::quantum::Pulse>("fc", channelName);
            xacc::InstructionParameter fcPreParameter(fcAngle);
            fcPre->setParameter(0, fcPreParameter);
            fcPre->setBits({0});

            auto fcPost = std::make_shared<xacc::quantum::Pulse>("fc", channelName);
            xacc::InstructionParameter fcPostParameter(-fcAngle);
            fcPost->setParameter(0, fcPostParameter);
            fcPost->setBits({0});
            fcPost->setStart(nSamples);

            cmddef_x_0->addInstructions({ fcPre, pulse, fcPost });
            xacc::contributeService(cmddef_x_0->name(), cmddef_x_0);
        }
       

        // ===============================================     

        
        // PI/2 rotation
        const int PULSE_PI_2_DURATION = 50;
        channelConfigs.addOrReplacePulse("X0_PI_2", QuaC::SquarePulse(PULSE_PI_2_DURATION));
        channelConfigs.addOrReplacePulse("X1_PI_2", QuaC::SquarePulse(PULSE_PI_2_DURATION));        
        setChannelConfigs(channelConfigs);

        // We will use the following pulse decomposition:
        // u3(θ, φ, λ) := U(θ, φ, λ) = Rz(φ + 3π)Rx(π/2)Rz(θ + π)Rx(π/2)Rz(λ)
        // i.e. two X(pi/2) pulses and a couple of frame changes

        // Add command defs:
        // Needs U3_0, U3_1, CX_0_1 and CX_1_0 
        {
            auto provider = xacc::getIRProvider("quantum");
            auto cmddef_u3_0 = provider->createComposite("pulse::u3_0");
            cmddef_u3_0->addVariables({"P0", "P1", "P2"});
            cmddef_u3_0->setBits({0});
            {
                auto fc_P0 = std::make_shared<xacc::quantum::Pulse>("fc", "d0");
                xacc::InstructionParameter fcParameterP0("P0 + pi");
                fc_P0->setParameter(0, fcParameterP0);
                fc_P0->setStart(PULSE_PI_2_DURATION);

                auto fc_P1 = std::make_shared<xacc::quantum::Pulse>("fc", "d0");
                xacc::InstructionParameter fcParameterP1("P1 + 3*pi");
                fc_P1->setParameter(0, fcParameterP1);
                fc_P1->setStart(PULSE_PI_2_DURATION + PULSE_PI_2_DURATION);

                auto fc_P2 = std::make_shared<xacc::quantum::Pulse>("fc", "d0");
                xacc::InstructionParameter fcParameterP2("P2");
                fc_P2->setParameter(0, fcParameterP2);
                fc_P2->setStart(0);
                
                auto firstPulseX_Pi_2 = std::make_shared<xacc::quantum::Pulse>("X0_PI_2", "d0");
                firstPulseX_Pi_2->setStart(0);
                
                auto secondPulseX_Pi_2 = std::make_shared<xacc::quantum::Pulse>("X0_PI_2", "d0");
                secondPulseX_Pi_2->setStart(PULSE_PI_2_DURATION);
                // Add pulse sequence to the cmd-def
                cmddef_u3_0->addInstructions({fc_P2, firstPulseX_Pi_2, fc_P0, secondPulseX_Pi_2, fc_P1});
            }

            auto cmddef_u3_1 = provider->createComposite("pulse::u3_1");
            cmddef_u3_1->addVariables({"P0", "P1", "P2"});
            cmddef_u3_1->setBits({1});
            {
                auto fc_P0 = std::make_shared<xacc::quantum::Pulse>("fc", "d1");
                xacc::InstructionParameter fcParameterP0("P0 + pi");
                fc_P0->setParameter(0, fcParameterP0);
                fc_P0->setStart(PULSE_PI_2_DURATION);

                auto fc_P1 = std::make_shared<xacc::quantum::Pulse>("fc", "d1");
                xacc::InstructionParameter fcParameterP1("P1 + 3*pi");
                fc_P1->setParameter(0, fcParameterP1);
                fc_P1->setStart(PULSE_PI_2_DURATION + PULSE_PI_2_DURATION);

                auto fc_P2 = std::make_shared<xacc::quantum::Pulse>("fc", "d1");
                xacc::InstructionParameter fcParameterP2("P2");
                fc_P2->setParameter(0, fcParameterP2);
                fc_P2->setStart(0);
                
                auto firstPulseX_Pi_2 = std::make_shared<xacc::quantum::Pulse>("X1_PI_2", "d1");
                firstPulseX_Pi_2->setStart(0);
                
                auto secondPulseX_Pi_2 = std::make_shared<xacc::quantum::Pulse>("X1_PI_2", "d1");
                secondPulseX_Pi_2->setStart(PULSE_PI_2_DURATION);
                // Add pulse sequence to the cmd-def
                cmddef_u3_1->addInstructions({fc_P2, firstPulseX_Pi_2, fc_P0, secondPulseX_Pi_2, fc_P1});
            }

            auto cmddef_cx_0_1 = provider->createComposite("pulse::cx_0_1");
            cmddef_cx_0_1->setBits({0, 1});
            
            auto cmddef_cx_1_0 = provider->createComposite("pulse::cx_1_0");
            cmddef_cx_1_0->setBits({1, 0});

            // TODO: define pulse libs and construct pulse sequence for the above 4 gates.  
            xacc::contributeService(cmddef_u3_0->name(), cmddef_u3_0);
            xacc::contributeService(cmddef_u3_1->name(), cmddef_u3_1);
            xacc::contributeService(cmddef_cx_0_1->name(), cmddef_cx_0_1);
            xacc::contributeService(cmddef_cx_1_0->name(), cmddef_cx_0_1);
        } 
    }

    const std::string name() const override { return "Default2Q"; }
};
}