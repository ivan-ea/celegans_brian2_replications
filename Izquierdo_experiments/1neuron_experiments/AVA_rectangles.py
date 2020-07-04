# -*- coding: utf-8 -*-
"""
@author: Ivan
"""

from brian2 import * 
from aux_functions import trimmed_triangle, trimmed_rectangle, save_plot, trimmed_rectangle_tail
import sys

n = 1
time_step = 0.0025*second
defaultclock.dt = time_step;
duration = 8000*ms

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
neurons.v = 0
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

store()

Imax = [30,20,10,5,2,-2,-5,-10,-20,-30]
trazas_I = []
trazas_V = []


for I in Imax:
    I_recorded = TimedArray((trimmed_rectangle_tail(duration*0.1/ms, I, duration)), dt=defaultclock.dt)
    restore()
    run(duration)
    trazas_V.append(trace.v[0])
    trazas_I.append(trace.I_i[0])
    
    subplot(211)
    plot(trace.t/ms, trace.v[0], label = 'v_'+str(I))
    
    subplot(212)
    plot(trace.t/ms, trace.I_i[0], label = 'I_'+str(I))

t_recorded = arange(int(duration/defaultclock.dt))*defaultclock.dt
xlabel('Time (ms)')
ylabel('i')
#legend(loc='upper right');
show()
sys.exit()


#plot(t_recorded, trace.v[0], label = 'v0')

save_plot("AVA_rectangle_"+str(duration),"t (ms) \t V_i \t I_i",t_recorded/ms,
          trazas_V[0],trazas_I[0],trazas_V[1],trazas_I[1],trazas_V[2],trazas_I[2],
          trazas_V[3],trazas_I[3],trazas_V[4],trazas_I[4],trazas_V[5],trazas_I[5],
          trazas_V[6],trazas_I[6],trazas_V[7],trazas_I[7],trazas_V[8],trazas_I[8],
          trazas_V[9],trazas_I[9])


sys.exit()
