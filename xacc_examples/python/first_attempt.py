# One-qubit pulse simulation: Rabi oscillation
# We need to have the XACC install directory in the Python path.
# Just in case users haven't already done that, set it here.
import sys
from pathlib import Path
sys.path.insert(1, str(Path.home()) + '/.xacc')

# Import xacc and quaC python wrapper
import os
import xacc
import json
import spectrum

import numpy as np
import matplotlib
# CADES VM don't have display
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# The Hamiltonian JSON object (OpenPulse format)
# omega0 = 2*pi, rotation speed: 100ns -> pi pulse (assume dt = 1) 
hamiltonianJson = {
        "description": "Hamiltonian of a one-qubit system.\n",
        "h_str": ["-0.5*omega0*Z0", "omegaa*X0||D0"],
        "osc": {},
        "qub": {
            "0": 2
        },
        "vars": {
            "omega0": 6.2831853,
            "omegaa": 0.0314159
        } 
}

# Create a pulse system model object 
model = xacc.createPulseModel()

# Load the Hamiltonian JSON (string) to the system model
loadResult = model.loadHamiltonianJson(json.dumps(hamiltonianJson))

qpu = xacc.getAccelerator('QuaC', {'system-model': model.name(), 'shots': 1024 })
channelConfig = xacc.BackendChannelConfigs()

T = 100
# Number of sample points to realize a PI pulse
nbSamples = 1000 #512 # 100
# dt (time between data samples)
channelConfig.dt = nbSamples / T #1.0
# omega0 = 2*pi => freq = 1.0 (drive at resonance) 
channelConfig.loFregs_dChannels = [1.0]
model.setChannelConfigs(channelConfig)

W = 0.02
k = int(2 * nbSamples * W)
n_orders = 36
# Initialize Slepians
Slepians, eigenvalues = spectrum.dpss(nbSamples, (nbSamples*W), k)
Slepians = Slepians[:, 0:n_orders]

state = 5 * np.ones(Slepians.shape[1])

# Square pulse with nbSamples elements
pulseData = (state * Slepians).sum(axis=1)
# Add that square pulse instruction to XACC
pulseName = 'Slepian'
xacc.addPulse(pulseName, pulseData)   
q = xacc.qalloc(1)
# Create the quantum program that contains the square pulse
# and the drive channel (D0) is set on the instruction
provider = xacc.getIRProvider('quantum')
prog = provider.createComposite('pulse')
slepianPulse = provider.createInstruction(pulseName, [0])
slepianPulse.setChannel('d0')
prog.addInstruction(slepianPulse)
# Measure Q0 (using the number of shots that was specified above)
prog.addInstruction(xacc.gate.create("Measure", [0]))

# Run the simulation
qpu.execute(q, prog)
resultProb = q.computeMeasurementProbability('1')

print(resultProb)

# Plot the result
plt.plot(pulseData)
os.chdir(os.path.dirname(os.path.abspath(__file__)))
plt.savefig('Slepian.png')