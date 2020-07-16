import numpy as np
import scipy.linalg as sp

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import sys, json
from pathlib import Path
sys.path.insert(1, str(Path.home()) + '/.xacc')
import xacc

nbSamples = 509
in_bW = 0.025
in_K = int(5) #int(2 * nbSamples * in_bW)
alpha_vector = np.ones(in_K, dtype=np.float64) 
x = np.array(xacc.SlepianPulse(alpha_vector, nbSamples, in_bW, in_K))

## DELETE FROM HERE
T = 100

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

model = xacc.createPulseModel()
# Load the Hamiltonian JSON (string) to the system model
loadResult = model.loadHamiltonianJson(json.dumps(hamiltonianJson))
qpu = xacc.getAccelerator('QuaC', {'system-model': model.name(), 'shots': 1024})
channelConfig = xacc.BackendChannelConfigs()
# Setting resolution of pulse
channelConfig.dt = nbSamples / T 
# Driving on resonance with qubit
channelConfig.loFregs_dChannels = [1.0]
model.setChannelConfigs(channelConfig)
pulseName = 'Slepian' 
xacc.addPulse(pulseName, x)
q = xacc.qalloc(1)
provider = xacc.getIRProvider('quantum')
prog = provider.createComposite('pulse')
slepianPulse = provider.createInstruction(pulseName, [0])
slepianPulse.setChannel('d0')
prog.addInstruction(slepianPulse)
prog.addInstruction(xacc.gate.create("Measure", [0]))
qpu.execute(q, prog)
resultProb = q.computeMeasurementProbability('1')
print(resultProb)

## TO HERE

plt.plot(x)
plt.savefig('first_order.png')