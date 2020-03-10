import os, json, sys, numpy as np
from pathlib import Path
sys.path.insert(1, str(Path.home()) + '/.xacc')
import xacc

# The Hamiltonian JSON object (OpenPulse format)
hamiltonianJson = {
    "description": "Hamiltonian of a one-qubit system.\n",
    "h_str": ["0.5*epsilon*Z0", "0.5*delta*X0||D0"],
    "osc": {},
    "qub": {
        "0": 2
    },
    "vars": {
        "epsilon": 0.0,
        "delta": 6.2831853
    }
}

# Create a pulse system model object
model = xacc.createPulseModel()

# Load the Hamiltonian JSON (string) to the system model
loadResult = model.loadHamiltonianJson(json.dumps(hamiltonianJson))

if loadResult is True:
    NB_SHOTS = 10000
    nSamples = 34.
    controlPulse = '0.062831853*exp(-t^2/(2*sigma^2))'
    # Run Quantum Optimal Control GOAT method to 
    # find gaussian pulse that approximates the X
    # gate on qubit 0
    goat = xacc.getOptimizer('quantum-control', 
                {'method':'GOAT', 'dimension':1, 
                    'target-U':'X0',
                    'control-params':['sigma'],
                    'control-funcs':[controlPulse],
                    'control-H':['X0'],
                    'max-time':nSamples,
                    'initial-parameters':[8.0]
                }
            )
    optimal_sigma = goat.optimize()[1][0]

    # Create the pulse at the optimal guassian sigma
    pulse = xacc.PulseFunc(controlPulse.replace('sigma', str(optimal_sigma)), int(nSamples))

    # Create the backend channel configuration
    channelConfigs = xacc.BackendChannelConfigs()
    channelConfigs.dt = 1
    channelConfigs.loFregs_dChannels = [0.0]
    channelConfigs.addOrReplacePulse('gaussian-x-gate', pulse)
    model.setChannelConfigs(channelConfigs)

    # Create the Pulse Program using the IRProvider
    provider = xacc.getIRProvider('quantum')
    compositeInst = provider.createComposite('pulse_qpt')
    pulseInst = xacc.createPulse('gaussian-x-gate', 'd0')
    pulseInst.setBits([0])
    compositeInst.addInstructions([pulseInst, provider.createInstruction('Measure',[0])])

    # Get reference to the QuaC Pulse simulator
    quaC = xacc.getAccelerator(
        'QuaC', {'system-model': model.name(), 'shots': NB_SHOTS})
     
    # Create the Quantum Process Tomography Algorithm
    qpt = xacc.getAlgorithm('qpt', {'circuit': compositeInst, 'accelerator': quaC, 'optimize-circuit': False})

    # Allocate a qubit and execute
    qubitReg = xacc.qalloc(1)
    qpt.execute(qubitReg)

    # Compute the fidelity with respect to the 
    # true process matrix for X
    chi_real_vec = [0., 0., 0., 0., 
                        0., 2., 0., 0., 
                        0., 0., 0., 0.,
                        0., 0., 0., 0.]
    result = qpt.calculate('fidelity', qubitReg, {
                           'chi-theoretical-real': chi_real_vec})
                        
    print('Fidelity = ', result)
    
    # Extra
    print(np.reshape(qubitReg['chi-real'],(4,4)))
    q = xacc.qalloc(1)
    quaC.execute(q, compositeInst)
    print(q)

    