#include "xacc.hpp"
#include <cassert>
#include "Pulse.hpp"
#include "PulseSystemModel.hpp"
#include "CommonGates.hpp"
#include "xacc_service.hpp"
#include <chrono>

int main (int argc, char** argv) {
    // Initialize the XACC Framework
    xacc::Initialize(argc, argv);
    // Load default pulse library for the IBM    
    const std::string jsonConfigFile = std::string(GATEIR_TEST_FILE_DIR) + "/paris_backend_pulse_lib.json";
    std::ifstream backendFile(jsonConfigFile);
    std::string jjson((std::istreambuf_iterator<char>(backendFile)), std::istreambuf_iterator<char>());
    
    // We cannot use the full Hamiltonian of the Paris device
    // because there are too many qubits (27) which requires an 
    // enormous amount of memory (esp. the superoperators)
    // Start with 2-q first, q0 and q1
    const std::string hamiltonianJson = R"#(
    {
        "description": "Qubits are modelled as a two level system.\n",
        "h_str": [
            "_SUM[i,0,1,wq{i}/2*(I{i}-Z{i})]",
            "_SUM[i,0,1,2*omegad{i}*X{i}||D{i}]",
            "jq0q1*Sp0*Sm1",
            "jq0q1*Sm0*Sp1",
            "2*omegad1*X0||U0",
            "2*omegad0*X1||U1"
        ],
        "osc": {},
        "qub": {
            "0": 2,
            "1": 2
        },
        "vars": {
            "jq0q1": 0.009864601,
            "jq10q12": 0.009990265,
            "jq11q14": 0.008859291,
            "jq12q13": 0.009864601,
            "jq12q15": 0.009676105,
            "jq13q14": 0.010053096,
            "jq14q16": 0.009676105,
            "jq15q18": 0.009047787,
            "jq16q19": 0.009236282,
            "jq17q18": 0.009424778,
            "jq18q21": 0.00948761,
            "jq19q20": 0.009424778,
            "jq19q22": 0.009299114,
            "jq1q2": 0.009864601,
            "jq1q4": 0.009236282,
            "jq21q23": 0.009676105,
            "jq22q25": 0.008796459,
            "jq23q24": 0.009047787,
            "jq24q25": 0.009047787,
            "jq25q26": 0.008984955,
            "jq2q3": 0.00678584,
            "jq3q5": 0.007979645,
            "jq4q7": 0.010555751,
            "jq5q8": 0.009361946,
            "jq6q7": 0.011309734,
            "jq7q10": 0.010115928,
            "jq8q11": 0.010304424,
            "jq8q9": 0.010555751,
            "omegad0": 0.41569,
            "omegad1": 0.3849,
            "omegad2": 0.3666,
            "omegad3": 0.45029,
            "omegad4": 0.38595,
            "omegad5": 0.34452,
            "omegad6": 0.25834,
            "omegad7": 0.31541,
            "omegad8": 0.33585,
            "omegad9": 0.35104,
            "omegad10": 0.32047,
            "omegad11": 0.28493,
            "omegad12": 0.40479,
            "omegad13": 0.4044,
            "omegad14": 0.39942,
            "omegad15": 0.34466,
            "omegad16": 0.48499,
            "omegad17": 0.45577,
            "omegad18": 0.44621,
            "omegad19": 0.34467,
            "omegad20": 0.366,
            "omegad21": 0.36219,
            "omegad22": 0.34883,
            "omegad23": 0.34314,
            "omegad24": 0.32014,
            "omegad25": 0.33578,
            "omegad26": 0.37613,
            "wq0": 31.869479461831972,
            "wq1": 31.544651525628915,
            "wq10": 30.91216305208345,
            "wq11": 31.157590888107645,
            "wq12": 31.682835699306896,
            "wq13": 32.13199165488961,
            "wq14": 30.777270490625337,
            "wq15": 30.464538672515012,
            "wq16": 31.54793041161031,
            "wq17": 31.73395864038912,
            "wq18": 31.017948154055063,
            "wq19": 29.92863798691419,
            "wq2": 30.275586236603104,
            "wq20": 31.655998333618296,
            "wq21": 30.42327453438617,
            "wq22": 31.357411955603823,
            "wq23": 32.206426833904096,
            "wq24": 31.284633452164687,
            "wq25": 30.490667365614655,
            "wq26": 31.15726353065678,
            "wq3": 30.756523103828073,
            "wq4": 31.954124895884757,
            "wq5": 30.142055653854847,
            "wq6": 32.68466503703604,
            "wq7": 32.30685481975422,
            "wq8": 31.847218476100185,
            "wq9": 32.481872017380006
        }
    })#";
    
    
    const std::vector<double> d_loFreqs {
        5.072185190116195,
        5.0204872184151395
    };
    
    const std::vector<double> u_loFreqs {
        5.0204872184151395,
        5.072185190116195
    };
    
    
    auto systemModel = std::make_shared<QuaC::PulseSystemModel>();
    const bool loadOk = systemModel->loadHamiltonianJson(hamiltonianJson);
    assert(loadOk);

    BackendChannelConfigs channelConfigs;
    channelConfigs.dt = 0.22222222222222221;
    channelConfigs.loFregs_dChannels = d_loFreqs;
    channelConfigs.loFregs_uChannels = u_loFreqs;

    systemModel->setChannelConfigs(channelConfigs);    
    // systemModel->setQubitInitialPopulation(0, 1.0);
    auto quaC = xacc::getAccelerator("QuaC", { std::make_pair("system-model", systemModel) });     
    // Contribute instruction
    quaC->contributeInstructions(jjson);     
    
    // TODO: the default pulse lib (from IBM) works *well* for single-qubit gates.
    // Two-qubit gates (i.e. CNOT) are not working as expected.
    // (it could be that the truncated model affects the cross-resonace dynamics
    // much more than single-qubit dynamics hence we don't get the expected CNOT gates
    // when using a sub-system model)
    auto compiler = xacc::getCompiler("xasm");
    auto ir = compiler->compile(R"(__qpu__ void f(qbit q) {
        H(q[0]);
        Z(q[0]);
        H(q[0]);
    })", nullptr);
    auto program = ir->getComposite("f");
    
    auto qubitReg = xacc::qalloc(2);
    quaC->execute(qubitReg, program);
    
    qubitReg->print();
    std::cout << "\n";
   
    // Finalize the XACC Framework
    xacc::Finalize();

    return 0;
}