import sys
from pathlib import Path
sys.path.insert(1, str(Path.home()) + '/.xacc')

# Import xacc and quaC python wrapper
import os
import xacc
import json
import numpy as np

# Three-level single-qubit Hamiltonian
# Ref: Phys. Rev. A 96, 022330 
# Equation (6), params in Sec III, first paragraph
hamiltonianJson = {
    "description": "One-qutrit Hamiltonian.",
    "h_latex": "",
    "h_str": ["(w - 0.5*alpha)*O0", "0.5*alpha*O0*O0", "O*(SM0 + SP0)||D0"],
    "osc": {},
    "qub": {
        "0": 3
    },
    "vars": {
        "w": 31.63772297724,
        "alpha": -1.47969,
        "O": 0.0314
    }
}

# Slepian Parameters
T = 100
nbSamples = 1000
in_bW = 0.02
in_K = int(2 * nbSamples * in_bW)
alpha_vector = np.ones(in_K)

# Create a pulse system model object 
model = xacc.createPulseModel()
# Load the Hamiltonian JSON (string) to the system model
loadResult = model.loadHamiltonianJson(json.dumps(hamiltonianJson))
qpu = xacc.getAccelerator('QuaC', {'system-model': model.name()})
channelConfig = xacc.BackendChannelConfigs()
# dt (time between data samples)
channelConfig.dt = nbSamples / T
# Drive at resonance: 31.63772297724/(2pi)    
channelConfig.loFregs_dChannels = [5.0353]

channelConfig.addOrReplacePulse('slepian', xacc.SlepianPulse(nbSamples, in_bW, in_K))
model.setChannelConfigs(channelConfig)

# Allocate qubits to run program
q = xacc.qalloc(1)
# Create the quantum program that contains the square pulse
# and the drive channel (D0) is set on the instruction
provider = xacc.getIRProvider('quantum')
prog = provider.createComposite('pulse')

# Calling xacc.createPulse() to create Slepian Matrix
slepian_matrix = xacc.createPulse('slepian_1', 'd0')
# Multiplying and summing by alpha vector to create pulse
pulse = (alpha_vector * slepian_matrix).sum(axis=1)
prog.addInstruction(pulse)

# Run the simulation
qpu.execute(q, prog)
# Get the probability of the |1> and |2> states
resultProb1[j][i] = q['DensityMatrixDiags'][1]
resultProb2[j][i] = q['DensityMatrixDiags'][2]
print(q)
