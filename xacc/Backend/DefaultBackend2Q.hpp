#pragma once

#include "xacc.hpp"
#include "xacc_service.hpp"
#include "PulseSystemModel.hpp"

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
                "h_str": ["0.5*(2*v0-alpha0)*O0", "0.5*alpha0*O0*O0", "r*X0||D0", "r*X0||U1", "r*X1||U0", "0.5*(2*v1-alpha1)*O1", "0.5*alpha1*O1*O1", "r*X1||D1", "j*(Sp0*Sm1+Sm0*Sp1)"],
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
                    "r": 0.125663706144,
                    "j": 0.062831853072
                }
            }
        )#";
        loadHamiltonianJson(hamiltonianJson);

        BackendChannelConfigs channelConfigs;
        channelConfigs.dt = 1.0;
        channelConfigs.loFregs_dChannels.emplace_back(4.9);
        channelConfigs.loFregs_dChannels.emplace_back(5.0);
        channelConfigs.loFregs_uChannels.emplace_back(4.9);
        channelConfigs.loFregs_uChannels.emplace_back(5.0);
        // Pulse library set-up:
        // TODO:
        channelConfigs.addOrReplacePulse("pulse1", {});
        setChannelConfigs(channelConfigs);

        // Add command defs:
        // Needs U3_0, U3_1, CX_0_1 and CX_1_0 
        {
            auto provider = xacc::getIRProvider("quantum");
            auto cmddef_u3_0 = provider->createComposite("pulse::u3_0");
            cmddef_u3_0->addVariables({"P0", "P1", "P2"});
            cmddef_u3_0->setBits({0});

            auto cmddef_u3_1 = provider->createComposite("pulse::u3_1");
            cmddef_u3_1->addVariables({"P0", "P1", "P2"});
            cmddef_u3_1->setBits({1});

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