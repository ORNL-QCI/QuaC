import sys, os, json, numpy as np
from pathlib import Path
sys.path.insert(1, str(Path.home()) + '/.xacc')
import xacc
import csv
import time
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

def runSim(qubitReg, controlState):
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
        
        # Set control qubit to 1 if set
        if controlState is True:
            model.setQubitInitialPopulation(0, 1.0)
       
        provider = xacc.getIRProvider('quantum')
        composite = provider.createComposite('test_pulse')
        pulse = xacc.createPulse('square', 'd0')
        composite.addInstruction(pulse)
        qpu.execute(qubitReg, composite)

# Run 1: control = 0    
# Two qutrits
qubitReg1 = xacc.qalloc(2)
runSim(qubitReg1, False)

# Wait for 1 second: in rare cases, the simulation runs so fast 
# that its output file timestamp (second resolution) may collide.
time.sleep(1)

# Run 2: control = 1    
# Two qutrits
qubitReg2 = xacc.qalloc(2)
runSim(qubitReg2, True)

# Get the CSV file name and data
csvFile1 = qubitReg1['csvFile']
data1 = np.genfromtxt(csvFile1, delimiter = ',', dtype=float, names=True)
csvFile2 = qubitReg2['csvFile']
data2 = np.genfromtxt(csvFile2, delimiter = ',', dtype=float, names=True)


# Replicate figs 3,5 in Sheldon PRA 2016
# Control qutrit: Pauli expectations for both control = 0 and 1
fig, ax = plt.subplots(6,1,sharey=True,figsize=(15,15))
ax[0].plot(data1['Time'], data1['X0'], 'b', label = '$|0\\rangle_C$')
ax[1].plot(data1['Time'], data1['Y0'], 'b', label = '$|0\\rangle_C$')
ax[2].plot(data1['Time'], data1['Z0'], 'b', label = '$|0\\rangle_C$')
ax[0].plot(data2['Time'], data2['X0'], 'r', label = '$|1\\rangle_C$')
ax[1].plot(data2['Time'], data2['Y0'], 'r', label = '$|1\\rangle_C$')
ax[2].plot(data2['Time'], data2['Z0'], 'r', label = '$|1\\rangle_C$')
ax[0].set_ylabel('$\\langle \hat{X}_C \\rangle$')
ax[1].set_ylabel('$\\langle \hat{Y}_C \\rangle$')
ax[2].set_ylabel('$\\langle \hat{Z}_C \\rangle$')

# Target qutrit: Pauli expectations for both control = 0 and 1
ax[3].plot(data1['Time'], data1['X1'], 'b', label = '$|0\\rangle_C$')
ax[4].plot(data1['Time'], data1['Y1'], 'b', label = '$|0\\rangle_C$')
ax[5].plot(data1['Time'], data1['Z1'], 'b', label = '$|0\\rangle_C$')
ax[3].plot(data2['Time'], data2['X1'], 'r', label = '$|1\\rangle_C$')
ax[4].plot(data2['Time'], data2['Y1'], 'r', label = '$|1\\rangle_C$')
ax[5].plot(data2['Time'], data2['Z1'], 'r', label = '$|1\\rangle_C$')
ax[3].set_ylabel('$\\langle \hat{X}_T \\rangle$')
ax[4].set_ylabel('$\\langle \hat{Y}_T \\rangle$')
ax[5].set_ylabel('$\\langle \hat{Z}_T \\rangle$')
ax[5].set_xlabel('Time ($10^{-10}s$)')
for i in range(6): 
    ax[i].set_xlim([0, 5100])
    ax[i].legend()

os.chdir(os.path.dirname(os.path.abspath(__file__)))
plt.savefig('RWA_pauli_dynamics.png')