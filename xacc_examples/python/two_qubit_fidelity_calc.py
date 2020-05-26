# Two-qubit default backend
import sys
from pathlib import Path
sys.path.insert(1, str(Path.home()) + '/.xacc')

# Import xacc and quaC python wrapper
import os, xacc, json
import numpy as np


# Get the two-qubit backend (QuaC default):
# This will load the default two-qubit backend:
# Backend Info (see DefaultBackend2Q.hpp for details):
"""
{
    "description": "Two-qubit Hamiltonian",
    "h_str": ["_SUM[i,0,1,wq{i}*O{i}]", "_SUM[i,0,1,delta{i}*O{i}*(O{i}-I{i})]", "_SUM[i,0,1,omegad{i}*X{i}||D{i}]", "omegad1*X0||U0", "omegad0*X1||U1", "jq0q1*Sp0*Sm1", "jq0q1*Sm0*Sp1"],
    "osc": {},
    "qub": {
        "0": 2,
        "1": 2
    },
    "vars": {
        "wq0": 30.518812656662774, 
        "wq1": 31.238229295532093,
        "delta0": -2.011875935,
        "delta1": -2.008734343,
        "omegad0": -1.703999855,
        "omegad1": -1.703999855,
        "jq0q1": 0.011749557 
    }
}
"""
# Format: <Accelerator>:<Backend name>
qpu = xacc.getAccelerator('QuaC:Default2Q')
# Example: drive a DRAG pulse on d0 channel which implements X gate on q0 
pulseName = "drag_pulse"
# Drag pulse sequence:
ampl = 0.165885
beta = 1.0
nSamples = 121
gamma = 20
pulseSequence = xacc.DragPulse(nSamples, ampl, gamma, beta)
#print("Pulse data: \n {}\n\n".format(pulseSequence))

# Register the pulse with XACC to construct pulse composite
# from these individual pulses:
xacc.addPulse(pulseName, pulseSequence)   

provider = xacc.getIRProvider('quantum')
prog = provider.createComposite('pulse_composite')

# Create pulse instruction and add to composite:
# We can add multiple pulse instructions to the composite.
dragPulseInst = provider.createInstruction(pulseName, [0])
dragPulseInst.setChannel('d0')
prog.addInstruction(dragPulseInst)

# 2-qubit register
qReg = xacc.qalloc(2)

# Density matrix fidelity calculation helper:
# The QuaC accelerator will calculate the fidelity if we provide
# a reference density matrix (expected density matrix)
# We need to provide the *real* and *imaginary* parts of the density matrix as real vectors. 
# Expected density matrix: rho = |10><10| since this drag pulse implements an X gate on the first qubit. 
expectedDmReal = np.array([
    0, 0, 0, 0,
    0, 0, 0, 0,
    0, 0, 1, 0,
    0, 0, 0, 0
], dtype = np.float64)
expectedDmImag = np.zeros(16)
          
# Add target density matrix info to the buffer before execution          
qReg.addExtraInfo("target-dm-real", expectedDmReal)
qReg.addExtraInfo("target-dm-imag", expectedDmImag)

# Run the simulation
qpu.execute(qReg, prog)
# print(qReg)

# Print the fidelity
fidelityResult = qReg["fidelity"]
print("\nFidelity: {}".format(fidelityResult))
