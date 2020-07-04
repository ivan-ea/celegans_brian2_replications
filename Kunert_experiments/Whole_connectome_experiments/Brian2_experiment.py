# -*- coding: utf-8 -*-
"""
Created on Sat May 30 17:39:38 2020

@author: Ivan Estevez Albuja
"""

from brian2 import *
from Wicks_model import Wicks_gap_jun, Wicks_chem_syn, Wicks_eqs_constIn, Wicks_eqs, calc_Vth
from Wicks_model import Ec, gs_u, gg_u, E_inh, s_eq, beta, t_syn, a_r, a_d, C, Gc, i_u

Syn_mat_original = loadtxt('connectome_data/Chem.csv', delimiter=',')
Gap_mat_original = loadtxt('connectome_data/Gap.csv', delimiter=',')
N = shape(Syn_mat_original)[0]

with open('connectome_data/labels_types.csv','r') as lf:
    labels_types = [line[0:-1].split(',') for line in lf] 
labels = [info[0] for info in labels_types]    
types =  [info[1] for info in labels_types] 

def idx_labels_start_with(start, labels):
    st_labels = list(filter(lambda x: x.startswith(start), labels))
    st_labels.sort()
    return [labels.index(l) for l in st_labels]

motor_indices   = [i for i, type_ in enumerate(types) if type_ == 'motor'] 
sensory_indices = [i for i, type_ in enumerate(types) if type_ == 'sensory']
inter_indices   = [i for i, type_ in enumerate(types) if type_ == 'inter']  

loco_colors = {'DA': 'orange', 'DB': 'red', 'DD': 'pink', 'VA': 'purple', 'VB': 'blue', 'VD': 'cyan','AS': 'green', }
loco_indices = {'DA': [], 'DB': [], 'DD': [], 'VA': [], 'VB': [], 'VD': []}
for label in loco_indices: 
    loco_indices[label] = idx_labels_start_with(label, labels) 
loco_indices['AS'] = [165,171,182,189,199,204,213,219,230,236,240]
loco_indices_rel = {}
for label in loco_indices: 
    loco_indices_rel[label] = [motor_indices.index(i) for i in loco_indices[label]]
all_loco_indices = concatenate(list(loco_indices.values()))

fwd_labels = ['DB', 'DD', 'VB', 'VD']
fwd_indices = {}
for l in fwd_labels: fwd_indices[l] = loco_indices[l].copy()
fwd_indices_rel = {}
for l in fwd_labels: fwd_indices_rel[l] = loco_indices_rel[l].copy()
all_fwd_indices = concatenate(list(fwd_indices.values()))    
   

# def ablation(adj_matrix, abl_indices):
#     return delete(delete(adj_matrix, abl_indices, axis = 0), abl_indices, axis = 1)
    
def brian2sim(t_sim,t_step,I_shape,rec_vars,in_indices,rec_indices,abl_indices):
    
    time_step = t_step*second
    defaultclock.dt = time_step
    duration = t_sim*second
    
    Syn_mat = loadtxt('connectome_data/Chem.csv', delimiter=',')
    Gap_mat = loadtxt('connectome_data/Gap.csv', delimiter=',')
        
    if abl_indices != None:
        for i in abl_indices:
            Syn_mat[i,:] = zeros(N)
            Syn_mat[:,i] = zeros(N)
            Gap_mat[i,:] = zeros(N)
            Gap_mat[:,i] = zeros(N)
            
    if type(I_shape) == int or type(I_shape) == float:
        is_I_const = True
        I_max = I_shape
        wicks_model = Wicks_eqs_constIn
    else:
        is_I_const = False
        I_max = I_shape[0]
        wicks_model = Wicks_eqs
        I_recorded = TimedArray(I_shape, dt=defaultclock.dt)
        
      
    neurons = NeuronGroup(N, model = wicks_model, method = 'rk2')
    neurons.v = Ec
    neurons.s = 0 
    neurons.v_th = Ec
    neurons.active_ = 0
    if(is_I_const): neurons.I_ext = 0
    for ix in in_indices: 
        neurons.active_[ix] = 1
        if(is_I_const): neurons.I_ext[ix] = I_max

      
        
    # Chemical Synapses    
    syn_sources, syn_targets = Syn_mat.nonzero()
    if len(syn_sources) > 0:
        S = Synapses(neurons, neurons, model = Wicks_chem_syn)
        S.connect(i = syn_sources, j = syn_targets)
        S.Gs = 0 * gs_u
        S.Esyn = 0 * mV
        
        for i, j in zip(syn_sources, syn_targets):
            weight = Syn_mat[i,j]
            S.Gs[i,j] = abs(weight) * gs_u
            if (weight < 0):
                S.Esyn[i,j] = E_inh
    
    
    # Gap Junctions   
    gap_sources, gap_targets = Gap_mat.nonzero()
    if len(gap_sources) > 0:
        G = Synapses(neurons, neurons, model = Wicks_gap_jun)
        G.connect(i = gap_sources, j = gap_targets)
        G.Gg = 0 * gg_u
        
        for i,j in zip(gap_sources, gap_targets):
            G.Gg[i,j] = Gap_mat[i,j] * gg_u
    
    
    # Calculation of thresholds
    thresholds = calc_Vth(N, Syn_mat, Gap_mat, I_max, in_indices)
    # print(thresholds)
    
    neurons.v_th = transpose(thresholds)
    
    trace = StateMonitor(source = neurons, variables=rec_vars, record = rec_indices)
    
    print('--Brian2 experiment starts:')
    print(N,'neurons,\t', t_sim, 'seconds,\t')
    print('Input neurons', [labels[i] for i in in_indices],' input current',str(I_max*i_u), '(constant: '+str(is_I_const)+')')
    if abl_indices != None: print('Ablated neurons', [labels[i] for i in abl_indices])
    print('Chemical synaptic connections: ', len(syn_sources), ', \t Gap junctions: ', len(gap_sources))
    print('time_step', t_step*1000, 'ms', ',\tvariables',rec_vars)
    print('Recording indices:',rec_indices)
    run(duration)
    print('--End of Brian2 experiment')
        
    return trace     