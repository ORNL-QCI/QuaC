#%%
import xacc
import gym
import spectrum
import json
#import gym_pulsecontrol
import numpy as np

from stable_baselines.common.policies import MlpPolicy
from stable_baselines import PPO2

# Total time, T, of control pulse
T = 100
# Number of pulse samples
nbSamples = 200
W = 0.05
k = int(2 * nbSamples * W)
n_orders = 15
# Initialize Slepians
Slepians, eigenvalues = spectrum.dpss(nbSamples, (nbSamples*W), k)
Slepians = Slepians[:, 0:n_orders]

gym.envs.register(
     id='PulseControl-v0',
     entry_point='gym_pulsecontrol.envs.pulsecontrol_env:PulseEnv',
)

env = gym.make('PulseControl-v0')
env.slepians_matrix = Slepians.copy()
env.n_orders = n_orders

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
env.model = xacc.createPulseModel()
# Load the Hamiltonian JSON (string) to the system model
loadResult = env.model.loadHamiltonianJson(json.dumps(hamiltonianJson))
env.qpu = xacc.getAccelerator('QuaC', {'system-model': env.model.name(), 'shots': 1024 })
env.channelConfig = xacc.BackendChannelConfigs()
# Setting resolution of pulse
env.channelConfig.dt = nbSamples / T 
# Driving on resonance with qubit
env.channelConfig.loFregs_dChannels = [1.0]
env.model.setChannelConfigs(env.channelConfig)

drl_model = PPO2('MlpPolicy', env,  
             verbose=0)
drl_model.learn(total_timesteps=50*n_orders)
