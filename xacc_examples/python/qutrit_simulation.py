# High-dimension system demonstration: qne-qutrit pulse simulation
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

# Create a pulse system model object 
model = xacc.createPulseModel()

# Load the Hamiltonian JSON (string) to the system model
loadResult = model.loadHamiltonianJson(json.dumps(hamiltonianJson))

if loadResult is True :
    qpu = xacc.getAccelerator('QuaC', {'system-model': model.name()})
    channelConfig = xacc.BackendChannelConfigs()
    # dt (time between data samples)
    channelConfig.dt = 1.0
    # Drive at resonance: 31.63772297724/(2pi)    
    channelConfig.loFregs_dChannels = [5.0353]

    # Number of sample points to realize a PI pulse
    nbSamples = 100
    
    # In the following example, we drive the system to a PI pulse
    # under various fictitious noise signal at 1->2 resonance.
    # Theoretical: the 1->2 transition is *detuned* from the
    # principle 0->1 transition by the aharmonicity (alpha = 1.47969)
    # Hence, any noise sources that have fourier component at that frequency
    # will cause leakage.
    # To demonstrate that, we use a square pulse + noise (@ alpha freq)
    # driving pulse (rather than a clean square pulse).
    # We observe the amount of leakage (prob. of |2> state) under various noise amplitude.

    # Vary the pulse width (0->pi rotation)   
    pulseWidth = np.linspace(1, nbSamples, 50)
    noiseAmp = np.array([0.0, 0.5, 1.0])
    resultProb1 = np.zeros((noiseAmp.size, pulseWidth.size))
    resultProb2 = np.zeros((noiseAmp.size, pulseWidth.size))

    j = 0
    for amp in noiseAmp:
        i = 0
        for width in pulseWidth:
            # Square pulse with nbSamples elements
            pulseData = np.ones(int(width))
            pulseName = 'pulseFn' + str(int(amp)) + str(int(width))

            # Add a square pulse + noise @ aharmonicity frequency
            # the noise has a configurable amplitude.
            # In practice, this noise term is decomposed (Fourier) from an arbitrary noise source.
            pulseFnStr = '1.0 + ' + str(amp) + '*cos(1.47969*t)'
            channelConfig.addOrReplacePulse(pulseName, xacc.PulseFunc(pulseFnStr, int(width), 1.0))
            model.setChannelConfigs(channelConfig)

            q = xacc.qalloc(1)
            # Create the quantum program that contains the square pulse
            # and the drive channel (D0) is set on the instruction
            provider = xacc.getIRProvider('quantum')
            prog = provider.createComposite('pulse')
            pulse = xacc.createPulse(pulseName, 'd0')
            prog.addInstruction(pulse)
            # Run the simulation
            qpu.execute(q, prog)
            # Get the probability of the |1> and |2> states
            resultProb1[j][i] = q['DensityMatrixDiags'][1]
            resultProb2[j][i] = q['DensityMatrixDiags'][2]
            print(q)
            i = i + 1
        j = j + 1
    
    # Plot the result
    fig, ax = plt.subplots(2, 1, sharex=True, figsize = (8, 5))
    plt.tight_layout()
    ax[0].plot(pulseWidth, resultProb1[0], 'b', label = '$Noise Amp = 0.0$')
    ax[0].plot(pulseWidth, resultProb1[1], 'g', label = '$Noise Amp = 0.5$')
    ax[0].plot(pulseWidth, resultProb1[2], 'r', label = '$Noise Amp = 1.0$')

    ax[1].plot(pulseWidth, resultProb2[0], 'b', label = '$Noise Amp = 0.0$')
    ax[1].plot(pulseWidth, resultProb2[1], 'g', label = '$Noise Amp = 0.5$')
    ax[1].plot(pulseWidth, resultProb2[2], 'r', label = '$Noise Amp = 1.0$')  

    ax[0].set_xlim([0, 100])
    ax[0].set_ylim([0.0, 1.0])
    ax[1].set_xlim([0, 100])
    ax[1].legend()
    ax[0].legend()
    ax[0].set_title('Probability of $|1\\rangle$ state')
    ax[1].set_title('Leakage - Probability of $|2\\rangle$ state')
    ax[0].set_ylabel('Prob($|1\\rangle$')
    ax[1].set_ylabel('Prob($|2\\rangle$')
    ax[1].set_xlabel('Time')
    plt.gcf().subplots_adjust(bottom=0.1)
    plt.gcf().subplots_adjust(left=0.1)
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    plt.savefig('Qutrit_Result.pdf')