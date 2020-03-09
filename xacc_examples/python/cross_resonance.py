import sys, os, json, numpy as np
from pathlib import Path
sys.path.insert(1, str(Path.home()) + '/.xacc')
import xacc
import csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# The Hamiltonian JSON object (OpenPulse format)
# Two-qutrit system in Target Resonator's Reference Frame (second qutrit) + RWA
# Reference: Sheldon, et. al., PRA 93, 060302(R) (2016)
hamiltonianJson = {
    "description": "Two-qutrit Hamiltonian. Rotating frame @ target qubit frame",
    "h_latex": "",
    "h_str": ["(w_0-w_1)*O0", "d*O0*(O0-I0)", "d*O1*(O1-I1)", "J*(SP0*SM1 + SM0*SP1)", "O*(SM0 + SP0)||D0"],
    "osc": {},
    "qub": {
        "0": 3,
        "1": 3
    },
    "vars": {
        "w_0": 5.114,
        "w_1": 4.914,
        "d": -0.33,
        "J": 0.004,
        "O": 0.060
    }
}

# Create a pulse system model object 
model = xacc.createPulseModel()

# Load the Hamiltonian JSON (string) to the system model
loadResult = model.loadHamiltonianJson(json.dumps(hamiltonianJson))

if loadResult is True :
    qpu = xacc.getAccelerator('QuaC', {'system-model': model.name(), 'shots': 1024, 'logging-period': 0.1 })
    
    nSamples = 1000
    channelConfigs = xacc.BackendChannelConfigs()
    channelConfigs.dt = 5.0
    # Since we are in the rotating frame (@ Q1 freq.), the LO freq of D0 (control qutrit)
    # is equivalent to zero (we drive Q0 @ the freq. of Q1)
    channelConfigs.loFregs_dChannels = [0.0]   
    channelConfigs.addOrReplacePulse('square', xacc.SquarePulse(nSamples))
    model.setChannelConfigs(channelConfigs)
    
    # Set control qubit to 1
    model.setQubitInitialPopulation(0, 1.0)
    # Two qutrits
    qubitReg = xacc.qalloc(2)

    provider = xacc.getIRProvider('quantum')
    composite = provider.createComposite('test_pulse')
    pulse = xacc.createPulse('square', 'd0')
    composite.addInstruction(pulse)
    qpu.execute(qubitReg, composite)
    
    # Get the CSV file name
    # TODO: reproduce Ipython notebook plots
    # This is just an example plot for demo atm.
    csvFile = qubitReg['csvFile']
    data = np.genfromtxt(csvFile, delimiter = ',', dtype=float, names=True)
    plt.scatter(data['Time'], data['Z1'])
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    plt.savefig('Plot.png')