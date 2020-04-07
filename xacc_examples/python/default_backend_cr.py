import sys
from pathlib import Path
sys.path.insert(1, str(Path.home()) + '/.xacc')

# This script is a Rabi-type simulation for two-qubit concurrence (cross-resonance)
# We vary the cross-resonance drive pulse width (square pulse) to achieve max entanglement.
# This is the pulse width to implement a [ZX]^1/2 for the CNOT.

# Import xacc and quaC python wrapper
import os
import xacc
import json
import numpy as np
import matplotlib
# CADES VM don't have display
matplotlib.use('Agg')
import matplotlib.pyplot as plt

qpu = xacc.getAccelerator('QuaC:Default2Q')

nbSamples = 100
pulseWidth = np.linspace(nbSamples, 10 * nbSamples, 10)
resultConcurrence = np.zeros(pulseWidth.size)

i = 0
for width in pulseWidth:
    # Square pulse with nbSamples elements
    pulseData = (1.0/ 0.0314) * np.ones(int(width))
    # Add that square pulse instruction to XACC
    pulseName = 'square' + str(width)
    xacc.addPulse(pulseName, pulseData)   
    q = xacc.qalloc(2)
    # Create the quantum program that contains the square pulse
    # and the U channel (CR drive) is set on the instruction
    provider = xacc.getIRProvider('quantum')
    prog = provider.createComposite('pulse')
    squarePulse = provider.createInstruction(pulseName, [0])
    squarePulse.setChannel('u1')
    prog.addInstruction(squarePulse)

    # Run the simulation
    concurrentInfoRequest = dict()
    concurrentInfoRequest['0_1'] = 0.0;
    q.addExtraInfo('Concurrence', concurrentInfoRequest)
    qpu.execute(q, prog)
    resultConcurrence[i] = q['Concurrence']['0_1']
    i = i + 1

#Plot the result
plt.scatter(pulseWidth, resultConcurrence)
os.chdir(os.path.dirname(os.path.abspath(__file__)))
plt.savefig('Concurrence_Rabi.png')
