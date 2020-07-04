# -*- coding: utf-8 -*-
"""
@author: Ivan Estevez Albuja
"""

from brian2 import * 
from Wicks_model import Wicks_gap_jun, Wicks_chem_syn, Wicks_eqs_constIn, calc_Vth
from aux_functions import constant_intput, diff_V, save_plot, visualise_connectivity
import sys

N = 5 
time_step = 0.0002*second
duration = 3*second
defaultclock.dt = time_step;


# Parameters (slight differences in Kunert and Jimin Kim)
C = 1.0*pF ;        Gc = 10*psiemens ;  Ec = -35.0*mV ; i_u = 1*pamp ;
gg_u = 100*psiemens ;  gs_u = 100*psiemens ;  E_inh = -48*mV ;
beta = 0.125/mV ;   t_syn = 1/second ;
a_r = 1.0/1.5 ;     a_d = 5.0/1.5 ;
#a_r = 1.0 ;     a_d = 5.0 ; #kunert
s_eq = a_r/(a_r+2*a_d)

# External input waveform
Imax = 1
I_shape = constant_intput(Imax,duration)

#Neurons of the Circuit
neurons = NeuronGroup(N, model =  Wicks_eqs_constIn, method = 'euler')
neurons.v = Ec 
neurons.s = 0 
neurons.v_th = 0
neurons.active_ = 0
ext_index = 0
neurons.active_[ext_index] = 1 # only neuron that will receive external input
neurons.I_ext = Imax #constant input

# Chemical Synapses
Syn_mat = array([[0, 1, -1, 1, 0],
                 [0, 0, 0, 0, 0],
                 [0, 0, 0, 0, 0],
                 [0, 0, 0, 0, 0],
                 [0, 0, 0, 0, 0]])

syn_sources, syn_targets = Syn_mat.nonzero()

S = Synapses(neurons, neurons, model = Wicks_chem_syn)
S.connect(i = syn_sources, j = syn_targets)
S.Gs = 0 * gs_u
S.Esyn = 0*mV

for i, j in zip(syn_sources, syn_targets):
    weight = Syn_mat[i,j]
    S.Gs[i,j] = abs(weight) * gs_u
    if (weight < 0):
        S.Esyn[i,j] = E_inh

#visualise_connectivity(S)

# Gap Junctions
Gap_mat = array([[0, 0, 0, 0, 0],
                 [0, 0, 0, 0, 0],
                 [0, 0, 0, 0, 0],
                 [0, 0, 0, 0, 1],
                 [0, 0, 0, 1, 0]])

gap_sources, gap_targets = Gap_mat.nonzero()

G = Synapses(neurons, neurons, model = Wicks_gap_jun)
G.connect(i = gap_sources, j = gap_targets)
G.Gg = 0 * gg_u

for i,j in zip(gap_sources, gap_targets):
    G.Gg[i,j] = Gap_mat[i,j] * gg_u

#visualise_connectivity(G)


neurons.v_th = transpose(calc_Vth(N,Syn_mat,Gap_mat,Gc,gg_u,gs_u,Ec,E_inh,s_eq,Imax,i_u,ext_index))



print(neurons.v_th)
#sys.exit()
#Input without recording I_ext variable
plot(arange(int(duration/defaultclock.dt))*defaultclock.dt, I_shape)
xlabel('Time (s)')
ylabel('I (pA)')
show()

trace = StateMonitor(neurons, ['v','s'], record=[0, 1, 2 ,3, 4])

run(duration)

#Plot results

for i in range(N):
    plot(trace.t/second, 1000*diff_V(trace[i].v, 0), label = 'V'+str(i))
    
xlabel('Time (s)')
ylabel('V (mV)')
legend(loc='upper right')
show()
for i in range(N):
    plot((trace.t/second)[4:], trace[i].s[4:], label = 's'+str(i))
    
xlabel('Time (s)')
ylabel('s')
legend(loc='upper right')
show()

sys.exit()
name = "K_5neuron_iconst_"+str(time_step/ms)+"_I"+str(Imax)+"_dur"+str(int(duration/second))
head="t (s) \t I (pA) \t V0 (mV) \t V1 (mV) \t V2 (mV) \t V3 (mV) \t V4 (mV)" 
head += " \t s0 \t s1 \t s2 \t s3 \t s4"
save_plot(name, head, trace.t, I_shape, 1000*diff_V(trace[0].v, 0),
          1000*diff_V(trace[1].v, 0), 1000*diff_V(trace[2].v, 0),
          1000*diff_V(trace[3].v, 0), 1000*diff_V(trace[4].v, 0),
          trace[0].s, trace[1].s, trace[2].s, trace[3].s, trace[4].s )
print("Written data in: ", name)
