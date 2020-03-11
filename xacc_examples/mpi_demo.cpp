#include "xacc.hpp"
#include <cassert>
#include "Pulse.hpp"
#include "PulseSystemModel.hpp"
#include "CommonGates.hpp"
#include "xacc_service.hpp"
#include <chrono>

// In this demo, we want to examine the MPI speed-up by varying the
// system dimension (2, 4, 6, 8, 10 qubits) and the number of MPI processes.
// In each test, pulses (X pulses) are applied to all qubits.

// Input params:
// -nq: number of qubits
// -nprocs: number of MPI processes
// This needs to be run by `mpiexec.hydra -n 1`, i.e. using MPI runtime but only request 1 process.
// additional processes are spawned internally.
// Run by: mpiexec.hydra -n 1 ./mpi_demo -nq <nQubits> -nprocs <nProcs>
int main (int argc, char** argv) {
    if (argc != 5)
    {      
        puts("Usage:\n \t-nq: Number of qubits (1 to 10).\n\t-nprocs: Number of MPI processes.");        
        return -1;
    }

    int nbQubits = -1;
    int nbProcs = -1;

    for(int i = 1; i < argc - 1; ++i)
    {
        if (std::string(argv[i]) == "-nq")
        {
            nbQubits = std::stoi(argv[i + 1]);
        }

        if (std::string(argv[i]) == "-nprocs")
        {
            nbProcs = std::stoi(argv[i + 1]);
        }
    } 

    if (nbQubits < 1 || nbQubits > 10)
    {
        puts("Invalid number of qubits.");        
        return -1;
    }

    if (nbProcs < 1)
    {
        puts("Invalid number of MPI processes.");        
        return -1;
    }

    // Create a Hamiltonian Json for a specific number of qubits
    const auto createHamiltonianJson = [](int in_nbQubits) {
        // A template Json to create Hamiltonian for an arbitrary number of qubits
        // Replace {{nbQubits}} with a number in [1, 9] range
        const std::string hamiltonianJsonTmpl = R"#(
        {
            "description": "Qubits are modelled as a two level system.\n",
            "h_str": ["_SUM[i,0,{{nbQubits}},wq{i}/2*Z{i}]", "_SUM[i,0,{{nbQubits}},omegad{i}*X{i}||D{i}]"],
            "osc": {},
            "qub": {
                "0": 2,
                "1": 2,
                "2": 2,
                "3": 2,
                "4": 2,
                "5": 2,
                "6": 2,
                "7": 2,
                "8": 2,
                "9": 2
            },
            "vars": {
                "omegad0": 1.25,
                "omegad1": 0.97, 
                "omegad2": 2.5,
                "omegad3": 1.05,
                "omegad4": 0.845,
                "omegad5": 1.25,
                "omegad6": 1.25,
                "omegad7": 0.97,
                "omegad8": 1.33,
                "omegad9": 1.45,
                "wq0": 30.91270129264568,
                "wq1": 30.36010168900955,
                "wq2": 31.041771660759178,
                "wq3": 28.36701429077905,
                "wq4": 29.298278199939336,
                "wq5": 31.14757431485366,
                "wq6": 31.387741224914162,
                "wq7": 30.232349262897678,
                "wq8": 31.502130591468386,
                "wq9": 31.769632280927663
            }
        })#";
        
        const std::string tmplPlaceholder = "{{nbQubits}}";
        std::string hamiltonianJson = hamiltonianJsonTmpl;
        // Replace both occurrences
        hamiltonianJson.replace(hamiltonianJson.find(tmplPlaceholder), tmplPlaceholder.length(), std::to_string(in_nbQubits - 1));
        hamiltonianJson.replace(hamiltonianJson.find(tmplPlaceholder), tmplPlaceholder.length(), std::to_string(in_nbQubits - 1));
        return hamiltonianJson;
    };


    // Initialize the XACC Framework
	xacc::Initialize(argc, argv);

    std::vector<double> loFreqs {
        30.91270129264568,
        30.36010168900955,
        31.041771660759178,
        28.36701429077905,
        29.298278199939336,
        31.14757431485366,
        31.387741224914162,
        30.232349262897678,
        31.502130591468386,
        31.769632280927663
    };

    for (auto& freq : loFreqs)
    {
        freq = freq/(2*M_PI);
    }

    auto systemModel = std::make_shared<QuaC::PulseSystemModel>();
    const bool loadOk = systemModel->loadHamiltonianJson(createHamiltonianJson(nbQubits));
    assert(loadOk);

    BackendChannelConfigs channelConfigs;
    channelConfigs.dt = 3.5555555555555554;
    channelConfigs.loFregs_dChannels = loFreqs;
    systemModel->setChannelConfigs(channelConfigs);    

    const int NB_SHOTS = 10000;
    const std::string mpiExecutorMode = "MPI::" + std::to_string(nbProcs);
    auto quaC = xacc::getAccelerator("QuaC", { std::make_pair("system-model", systemModel), 
                                                std::make_pair("shots", NB_SHOTS),
                                                std::make_pair("execution-mode", mpiExecutorMode) });


    // Contribute all cmd-defs in the backend Json as XACC instruction.
    // This will activate Gates -> Pulses decomposition when simulating the circuit.
    const std::string jsonConfigFile = std::string(GATEIR_TEST_FILE_DIR) + "/test_backends.json";
    std::cout << "File name: " << jsonConfigFile << "\n";
    std::ifstream backendFile(jsonConfigFile);
    std::string jjson((std::istreambuf_iterator<char>(backendFile)), std::istreambuf_iterator<char>());
    quaC->contributeInstructions(jjson);     

    // Note: this may take a long time to run.
    auto qubitReg = xacc::qalloc(nbQubits);
    auto provider = xacc::getIRProvider("quantum");
    auto program = provider->createComposite("test_pulse");

    // List of Q0 -> Q9 pulses that we may use
    const std::vector<std::string> pulseList {
        "Xp_d0_ddeb",
        "Xp_d1_9f56",
        "Xp_d2_96dd",
        "Xp_d3_bb42",
        "Xp_d4_d499",
        "Xp_d5_adb5",
        "Xp_d6_93e1",
        "Xp_d7_e6f5",
        "Xp_d8_683a",
        "Xp_d9_e4e5"
    };

    for (int i = 0; i < nbQubits; ++i)
    {
        auto pulseInst = xacc::getContributedService<xacc::Instruction>(pulseList[i]);
        pulseInst->setChannel("d" + std::to_string(i));
        program->addInstruction(pulseInst);
    }

    for (int i = 0; i < nbQubits; ++i)
    {
        // Add a measure to get bitstrings
  	    auto meas = std::make_shared<xacc::quantum::Measure>(i);
        program->addInstruction(meas);
    }
    
    // Time the execution 
    const std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
   
    // Run the Pulse simulation with the Hamiltonian provided
    quaC->execute(qubitReg, program);
    const std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    
    std::cout << "Qubits = " << nbQubits << ", Procs = " << nbProcs << ": Total elapsed time = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << " [ms]\n";

    qubitReg->print();
    std::cout << "\n";

    // Finalize the XACC Framework
    xacc::Finalize();
    return 0;
}
