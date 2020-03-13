import sys, os, json, numpy as np
from pathlib import Path
sys.path.insert(1, str(Path.home()) + '/.xacc')
import xacc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# The Hamiltonian JSON object (OpenPulse format)
# omega0 = 2*pi, rotation speed: 100ns -> pi pulse (assume dt = 1) 
hamiltonianJson = {
        "description": "Hamiltonian of a one-qubit system.\n",
        "h_str": ["omega0*Z0", "omegaa*X0||D0"],
        "osc": {},
        "qub": {
            "0": 2
        },
        "vars": {
            "omega0": 0.00,
            "omegaa": 0.02
        } 
}

# Create a pulse system model object 
model = xacc.createPulseModel()

# Load the Hamiltonian JSON (string) to the system model
loadResult = model.loadHamiltonianJson(json.dumps(hamiltonianJson))

if loadResult is True :
    qpu = xacc.getAccelerator('QuaC', {'system-model': model.name(), 'shots': 1024, 'logging-period': 0.1 })
    
    nSamples = 1000
    channelConfigs = xacc.BackendChannelConfigs()
    channelConfigs.dt = 0.1
    channelConfigs.loFregs_dChannels = [1e-128]

    fourierSeries = '''
        0.726924 + 
        0.065903*cos(1*0.1*t) + 0.128627*sin(1*0.1*t) +
        0.079360*cos(2*0.1*t) + 0.111686*sin(2*0.1*t) + 
        0.096717*cos(3*0.1*t) + 0.096822*sin(3*0.1*t) + 
        0.106937*cos(4*0.1*t) + 0.092216*sin(4*0.1*t) + 
        0.215306*cos(5*0.1*t) + 0.118562*sin(5*0.1*t) +
        0.117682*cos(6*0.1*t) + 0.126134*sin(6*0.1*t) + 
        0.100447*cos(7*0.1*t) + 0.120409*sin(7*0.1*t) + 
        0.103292*cos(8*0.1*t) + 0.108712*sin(8*0.1*t)'''
    
    channelConfigs.addOrReplacePulse('fourier', xacc.PulseFunc(fourierSeries, nSamples, channelConfigs.dt))
    model.setChannelConfigs(channelConfigs)

    qubitReg = xacc.qalloc(1)

    provider = xacc.getIRProvider('quantum')
    composite = provider.createComposite('test_pulse')
    pulse = xacc.createPulse('fourier', 'd0')
    composite.addInstruction(pulse)

    qpu.execute(qubitReg, composite)
    print(qubitReg)
    # Retrieve time-stepping raw data
    csvFile = qubitReg['csvFile']
    data = np.genfromtxt(csvFile, delimiter = ',', dtype=float, names=True)
    fig, ax = plt.subplots(2, 1, sharex=True, figsize = (8, 5))
    plt.tight_layout()
    ax[0].plot(data['Time'], data['Channel0'], 'b', label = '$D_0(t)$')
    ax[1].plot(data['Time'], data['X0'], 'b', label = '$\\langle X \\rangle$')
    ax[1].plot(data['Time'], data['Y0'], 'g', label = '$\\langle Y \\rangle$')
    ax[1].plot(data['Time'], data['Z0'], 'r', label = '$\\langle Z \\rangle$')
    ax[1].plot(data['Time'], data['Population0'], 'k', label = '$Prob(1)$')

    ax[0].set_xlim([0, 100])
    ax[0].set_ylim([-0.1, 2.0])
    ax[1].set_xlim([0, 100])
    ax[1].set_ylim([-1.1, 1.1])
    ax[1].legend()
    ax[1].legend(loc='upper center', bbox_to_anchor=(0.5, 1.15), fancybox=True, shadow=True, ncol=5)
    ax[0].legend()
    ax[1].set_xlabel('Time')
    plt.gcf().subplots_adjust(bottom=0.1)
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    plt.savefig('Fourier_Pulse_Response.pdf')