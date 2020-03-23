import sys, os, json, numpy as np
from pathlib import Path
sys.path.insert(1, str(Path.home()) + '/.xacc')
import xacc

# The Hamiltonian JSON object (OpenPulse format)
hamiltonianJson = {
        "description": "Hamiltonian of a one-qubit system.\n",
        "h_str": ["-0.5*omega0*Z0", "omegaa*X0||D0"],
        "osc": {},
        "qub": {
            "0": 2
        },
        "vars": {
            "omega0": 0.0,
            "omegaa": 0.062832
        } 
}

# Create a pulse system model object 
model = xacc.createPulseModel()

# Load the Hamiltonian JSON (string) to the system model
loadResult = model.loadHamiltonianJson(json.dumps(hamiltonianJson))

if loadResult is True :
    qpu = xacc.getAccelerator('QuaC', {'system-model': model.name()})    
    channelConfigs = xacc.BackendChannelConfigs()
    channelConfigs.dt = 1.0
    channelConfigs.loFregs_dChannels = [0.0]
    model.setChannelConfigs(channelConfigs)

    # Get the XASM compiler
    xasmCompiler = xacc.getCompiler('xasm');
    # Composite to be transform to pulse
    ir = xasmCompiler.compile('''__qpu__ void f(qbit q) {
        Rx(q[0], 1.57);
    }''', qpu);
    program = ir.getComposites()[0]

    # Run the pulse IRTransformation 
    optimizer = xacc.getIRTransformation('quantum-control')
    optimizer.apply(program, qpu, {
        'method': 'GOAT',
        'control-params': ['sigma'],
        # Gaussian pulse
        'control-funcs': ['exp(-t^2/(2*sigma^2))'],
        # Initial params
        'initial-parameters': [8.0],
        'max-time': 100.0
    })

    # This composite should be a pulse composite now
    print(program)
    
    # Run the simulation of the optimized pulse program
    qubitReg = xacc.qalloc(1)
    qpu.execute(qubitReg, program)
    print(qubitReg)