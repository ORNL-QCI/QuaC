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

    # Get Quantum Process Tomography Algo
    qpt = xacc.getAlgorithm('qpt')

    # Compute Theoretical Chi Matrix
    qppReg = xacc.qalloc(1)
    acc = xacc.getAccelerator('qpp', {'shots': NB_SHOTS})
    compiler = xacc.getCompiler('xasm')
    ir = compiler.compile('''__qpu__ void f(qbit q) {
        X(q[0]);
    }''', None)

    qppCompositeInstr = ir.getComposites()[0]
    qpt.initialize({'circuit': qppCompositeInstr, 'accelerator': acc})
    qpt.execute(qppReg)
    qpp_chi_real_vec = qppReg["chi-real"]
    qpp_chi_imag_vec = qppReg["chi-imag"]

    # Setup Backend Channel Configuration
    channelConfigs = xacc.BackendChannelConfigs()
    channelConfigs.dt = 0.005
    channelConfigs.loFregs_dChannels = [0.0]
    nbSamples = 100
    channelConfigs.addOrReplacePulse('square', xacc.SquarePulse(nbSamples))
    model.setChannelConfigs(channelConfigs)

    # Create the Pulse Program, a Square Pulse
    qubitReg = xacc.qalloc(1)
    provider = xacc.getIRProvider('quantum')
    compositeInst = provider.createComposite('pulse_qpt')
    pulseInst = xacc.createPulse('square', 'd0')
    pulseInst.setBits([0])
    compositeInst.addInstruction(pulseInst)

    # Execute on QuaC Pulse Simulator
    quaC = xacc.getAccelerator(
        'QuaC', {'system-model': model.name(), 'shots': NB_SHOTS})
    qpt.initialize({'circuit': compositeInst, 'accelerator': quaC, 'optimize-circuit': False})
    qpt.execute(qubitReg)

    # Results stored on qubitReg, use 
    # it to compute the fidelity now
    
    # Compute the Fidelity
    result = qpt.calculate('fidelity', qubitReg, {
                           'chi-theoretical-real': qpp_chi_real_vec, 'chi-theoretical-imag': qpp_chi_imag_vec})
    print('Fidelity = ', result)
