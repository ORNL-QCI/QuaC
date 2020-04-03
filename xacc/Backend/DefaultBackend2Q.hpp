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
                "description": "Two-qubit Hamiltonian.",
                "h_latex": "",
                "h_str": ["0.5*(2*v0-alpha0)*O0", "0.5*alpha0*O0*O0", "r*(SM0 + SP0)||D0", "r*(SM0 + SP0)||U1", "r*(SM1 + SP1)||U0", "0.5*(2*v1-alpha1)*O1", "0.5*alpha1*O1*O1", "r*(SM1 + SP1)||D1", "j*(Sp0*Sm1+Sm0*Sp1)"],
                "osc": {},
                "qub": {
                    "0": 2,
                    "1": 2
                },
                "vars": {
                    "v0": 31.415926535898,
                    "v1": 32.044245066616,
                    "alpha0": -2.073451151369, 
                    "alpha1": -2.073451151369,
                    "r": 0.0314,
                    "j": 0.0062831853072
                }
            }
        )#";
        loadHamiltonianJson(hamiltonianJson);

        BackendChannelConfigs channelConfigs;
        channelConfigs.dt = 1.0;
        channelConfigs.loFregs_dChannels.emplace_back(5.0);
        channelConfigs.loFregs_dChannels.emplace_back(5.1);
        channelConfigs.loFregs_uChannels.emplace_back(5.0);
        channelConfigs.loFregs_uChannels.emplace_back(5.1);
        // Pulse library set-up:
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