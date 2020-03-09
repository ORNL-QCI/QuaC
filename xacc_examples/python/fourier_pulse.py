import sys, os, json, numpy as np
from pathlib import Path
sys.path.insert(1, str(Path.home()) + '/.xacc')
import xacc

# The Hamiltonian JSON object (OpenPulse format)
# omega0 = 2*pi, rotation speed: 100ns -> pi pulse (assume dt = 1) 
hamiltonianJson = {
        "description": "Hamiltonian of a one-qubit system.\n",
        "h_str": ["omega0*Z0", "omegaa*X0||D0"],
        "osc": {},
        "qub": {
            "0": 2
        },
        "vars": {
            "omega0": 0.0,
            "omegaa": 0.02
        } 
}

# Create a pulse system model object 
model = xacc.createPulseModel()

# Load the Hamiltonian JSON (string) to the system model
loadResult = model.loadHamiltonianJson(json.dumps(hamiltonianJson))

if loadResult is True :
    qpu = xacc.getAccelerator('QuaC', {'system-model': model.name(), 'shots': 1024 })
    
    nSamples = 1000
    channelConfigs = xacc.BackendChannelConfigs()
    channelConfigs.dt = 0.1
    channelConfigs.loFregs_dChannels = [0.0]

    fourierSeries = '''
        0.726924 + 
        0.065903*cos(1*0.1*t) + 0.128627*sin(1*0.1*t) +
        0.079360*cos(2*0.1*t) + 0.111686*sin(2*0.1*t) + 
        0.096717*cos(3*0.1*t) + 0.096822*sin(3*0.1*t) + 
        0.106937*cos(4*0.1*t) + 0.092216*sin(4*0.1*t) + 
        0.215306*cos(5*0.1*t) + 0.118562*sin(5*0.1*t) +
        0.117682*cos(6*0.1*t) + 0.126134*sin(6*0.1*t) + 
        0.100447*cos(7*0.1*t) + 0.120409*sin(7*0.1*t) + 
        0.103292*cos(8*0.1*t) + 0.108712*sin(8*0.1*t)'''
    
    channelConfigs.addOrReplacePulse('fourier', xacc.PulseFunc(fourierSeries, nSamples, channelConfigs.dt))
    model.setChannelConfigs(channelConfigs)

    qubitReg = xacc.qalloc(1)

    provider = xacc.getIRProvider('quantum')
    composite = provider.createComposite('test_pulse')
    pulse = xacc.createPulse('fourier', 'd0')
    composite.addInstruction(pulse)

    qpu.execute(qubitReg, composite)
    print(qubitReg)
