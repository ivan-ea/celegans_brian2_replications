# -*- coding: utf-8 -*-
"""
Created on Sun Apr 19 10:03:20 2020

@author: Ivan Estevez Albuja
"""

from brian2 import * 
from aux_functions import  save_plot, visualise_connectivity
import sys

n = 3

#tau = 10*ms 
time_step = 0.0025*second
defaultclock.dt = time_step;


eqsIzq = '''
dv/dt = (I_i - v + Igap + Isyn) / tau_i : 1
Igap : 1 # gap junction current
Isyn : 1 # chemical synapsis current
# per neuron parameters
I_i : 1 (constant)
tau_i : second (constant)
bias_i : 1 (constant)
label : integer (constant)
'''

gap_jun_eqs = '''
g : 1 # gap junction 'conductance', always > 0
Igap_post = g * (v_pre - v_post) : 1 (summed)
'''

syn_eqs = '''
Isyn_post = w/(1+exp(-(v_pre + bias_i_pre))): 1 (summed)
w : 1
'''

neurons = NeuronGroup(n, eqsIzq, #threshold = 'v > 1', reset = 'v = 0', #no th no reset => no spiking
                      method = 'rk4')

AS, DA, DB = 0, 1, 2
neurons.label = [AS, DA, DB]
#initial conditions for v, for every neuron
duration = 8000*ms
neurons.v = [-13,0,0]
neurons.I_i = 0 #[0, 0.3, 0.5]
neurons.tau_i = [0.05, 0.15, 0.2]*second
neurons.bias_i = [4.4, -1.8, -3.3]

#chemical synapses
S = Synapses(neurons, neurons, model=syn_eqs)
S.connect(i=[0,0,1,1,2,2], j=[0,1,1,2,2,0])
S.w[0,0] = '7.9'
S.w[0,1] = '2.5'
S.w[1,1] = '3'
S.w[1,2] = '2.5' #1.5 exponential, 2.5 periodic
S.w[2,2] = '3.0'
S.w[2,0] = '-18.3'

trace = StateMonitor(neurons, 'v', record=[0, 1, 2])


#gap junction connections
GapJ = Synapses(neurons, neurons, gap_jun_eqs)
# GapJ.connect()
# GapJ.g = .02

run(duration)

plot(trace.t/ms, trace[0].v, label = '0')
plot(trace.t/ms, trace[1].v, label = '1')
plot(trace.t/ms, trace[2].v, label = '2')
save_plot("3_neuron_"+"%0.4f"%(abs(duration/ms*mean(S.w))), 
          "t (ms) \t V_as \t V_da \t V_db",
          trace.t/ms,trace[0].v,trace[1].v,trace[2].v)
xlabel('Time (ms)')
ylabel('v')
legend(loc='upper right');
show()

sys.exit()


visualise_connectivity(S)