# -*- coding: utf-8 -*-
"""
@author: Ivan Estevez Albuja
"""

from brian2 import * 
from aux_functions import trimmed_triangle, trimmed_rectangle, save_plot
import sys

n = 1
time_step = 0.0025*second
defaultclock.dt = time_step;
duration = 4000*ms

Imax = 30

t_recorded = arange(int(duration/defaultclock.dt))*defaultclock.dt
#I_recorded = TimedArray((trimmed_rectangle(duration*0.1/ms, 20,duration)), dt=defaultclock.dt)
I_recorded = TimedArray((trimmed_triangle(duration*0.1/ms,duration, Imax)), dt=defaultclock.dt)


eqsIzq = '''
dv/dt = (active*I_i - v + Igap + Isyn) / tau_i : 1
Igap : 1 # gap junction current
Isyn : 1 # chemical synapsis current
# per neuron parameters
active : 1
I_i = I_recorded(t) : 1
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

neurons = NeuronGroup(n, model= eqsIzq, method = 'euler')


#initial conditions for v, for every neuron
neurons.v = -30
neurons.active_ = 0
neurons.active_[0] = 1
neurons.tau_i = 0.03*second
neurons.bias_i = -1.8

trace = StateMonitor(neurons, variables = True, record=0)

#chemical synapses
S = Synapses(neurons, neurons, model=syn_eqs)
S.connect(i=[0], j=[0])
S.w[0,0] = '0.6'

#gap junction connections
#GapJ = Synapses(neurons, neurons, gap_jun_eqs)
# GapJ.connect()
# GapJ.g = .02

run(duration)



plot(trace.t/ms, trace.v[0], label = 'v')
ylabel('v')
xlabel('Time (ms)')
legend(loc='upper right');
show()
plot(trace.t/ms, trace.I_i[0], label = 'I')
save_plot("AVA_triangle_"+str(duration),"t (ms) \t V \t I",trace.t/ms,trace.v[0],trace.I_i[0])
xlabel('Time (ms)')
ylabel('I')
legend(loc='upper right');
show()

sys.exit()
