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
# Note to ORNL CADES users:
# Make sure numpy and matplotlib are installed
# If not, these can be installed by:
# sudo apt install python3-numpy
# sudo apt install python3-matplotlib
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

if loadResult is True :
    qpu = xacc.getAccelerator('QuaC', {'system-model': model.name(), 'shots': 1024 })
    channelConfig = xacc.BackendChannelConfigs()
    # dt (time between data samples)
    channelConfig.dt = 1.0
    # omega0 = 2*pi => freq = 1.0 (drive at resonance) 
    channelConfig.loFregs_dChannels = [1.0]
    model.setChannelConfigs(channelConfig)

    # Number of sample points to realize a PI pulse
    nbSamples = 100
    
    # In the followings, we run a Rabi simulation by changing the pulse width
    # hence we should have a sinusoidal respond of the probability the output is measured to be 1.
    # Run up to 4*pi rotation (2 full Rabi cycle) 
    pulseWidth = np.linspace(1, 4 * nbSamples, 50)
    resultProb = np.zeros(pulseWidth.size)
    
    i = 0
    for width in pulseWidth:
        # Square pulse with nbSamples elements
        pulseData = np.ones(int(width)
        # Add that square pulse instruction to XACC
        pulseName = 'square' + str(width)
        xacc.addPulse(pulseName, pulseData)   
        q = xacc.qalloc(1)
        # Create the quantum program that contains the square pulse
        # and the drive channel (D0) is set on the instruction
        provider = xacc.getIRProvider('quantum')
        prog = provider.createComposite('pulse')
        squarePulse = provider.createInstruction(pulseName, [0])
        squarePulse.setChannel('d0')
        prog.addInstruction(squarePulse)
        # Measure Q0 (using the number of shots that was specified above)
        prog.addInstruction(xacc.gate.create("Measure", [0]))

        # Run the simulation
        qpu.execute(q, prog)
        resultProb[i] = q.computeMeasurementProbability('1')
        i = i + 1
    
    # Plot the result
    plt.scatter(pulseWidth, resultProb)
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    plt.savefig('Rabi.png')
else :
    print("Failed to load Hamiltonian Json. Please check your input.")
