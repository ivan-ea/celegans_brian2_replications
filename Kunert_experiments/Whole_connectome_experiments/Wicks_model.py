# -*- coding: utf-8 -*-
"""
Created on Mon May 11 11:31:45 2020

@author: Ivan Estevez Albuja

# Neuron model from C. elegans by Wicks et. al.
# Original publication: https://doi.org/10.1523/JNEUROSCI.16-12-04017.1996 
"""

from brian2 import * 


# Parameters (slight differences in Kunert and Jimin Kim papers)

C = 1.0*pF ;        Gc = 10*psiemens ;  Ec = -35.0*mV ; i_u = 1*pamp ;
gg_u = 100*psiemens ;  gs_u = 100*psiemens ;  E_inh = -45*mV ;

#a_r = 1.0/1.5 ;     a_d = 5.0/1.5 ; # jimin
a_r = 1.0 ;     a_d = 5.0 ; # kunert
beta = 0.125/mV ;   t_syn = 1/second ;
s_eq = a_r/(a_r+2*a_d)

# Equations 

Wicks_eqs_noInput = '''
dv/dt = (active*I_ext*i_u - Gc*(v - Ec) - Igap - Isyn) / C : volt
Igap : ampere # gap junction current
Isyn : ampere # chemical synapsis current
ds/dt = t_syn * a_r * (1 - s) / (1 + exp(-beta * (v - v_th))) - t_syn * a_d * s : 1
# per neuron parameters
v_th : volt
active : 1
'''

Wicks_eqs = Wicks_eqs_noInput + 'I_ext = I_recorded(t) : 1' + ' \n'

Wicks_eqs_constIn = Wicks_eqs_noInput + 'I_ext : 1 (constant)' + ' \n'

Wicks_gap_jun = '''
Gg : siemens # gap junction 'conductance', always > 0 (from the matrix)
Igap_post = Gg  * (v_post - v_pre) : ampere (summed)
'''

Wicks_chem_syn = '''
Gs : siemens # connection weight (from the matrix)
Esyn : volt #  0mV for excitatory synapses and −45 mV for inhibitory synapses
Isyn_post = Gs  * s_pre * (v_post - Esyn): ampere (summed)
'''

# calculate vthreshold

def calc_Vth(N, Chem_m, Gap_m, Iext, input_indices):  
   
    M1 = eye(N)*Gc
    
    # gap junctions contribution
    M2 = zeros([N,N])*gg_u
    gap_sources, gap_targets = Gap_m.nonzero()
    for i, j in zip(gap_sources, gap_targets):
        M2[i,j] = -Gap_m[i,j]*gg_u
        M2[i,i] += Gap_m[i,j]*gg_u
    
    # chemical synapses contribution
    M3 = zeros([N,N])*gs_u
    chem_sources, chem_targets = Chem_m.nonzero()
    for i, j in zip(chem_sources, chem_targets):
        M3[j,j] += s_eq*abs(Chem_m[i,j])*gs_u
    
    A = M1 + M2 + M3


    b1 = ones(N)*Gc*Ec
    
    # input current
    b2 = zeros(N)*i_u
    for ix in input_indices:
        b2[ix] = Iext*i_u
    
    # inhibitory synapses
    b3 = zeros(N)*i_u
    for i, j in zip(chem_sources, chem_targets):
        if (Chem_m[i,j] < 0):
            b3[j] += s_eq*abs(Chem_m[i,j])*gs_u  * E_inh

    b = b1+b2+b3
    

    # linear system Av = b
    v_th = linalg.solve(A, transpose(b)) * volt

    return v_th

def calc_Vth_1input(N, Chem_m, Gap_m, Iext, ext_index):     
    M1 = eye(N)*Gc
    
    # gap junctions contribution
    M2 = zeros([N,N])*gg_u
    gap_sources, gap_targets = Gap_m.nonzero()
    for i, j in zip(gap_sources, gap_targets):
        M2[i,j] = -Gap_m[i,j]*gg_u
        M2[i,i] += Gap_m[i,j]*gg_u
    
    # chemical synapses contribution
    M3 = zeros([N,N])*gs_u
    chem_sources, chem_targets = Chem_m.nonzero()
    for i, j in zip(chem_sources, chem_targets):
        M3[j,j] += s_eq*abs(Chem_m[i,j])*gs_u
    
    A = M1 + M2 + M3


    b1 = ones(N)*Gc*Ec
    
    # input current
    b2 = zeros(N)*i_u
    b2[ext_index] = Iext*i_u
    
    # inhibitory synapses
    b3 = zeros(N)*i_u
    for i, j in zip(chem_sources, chem_targets):
        if (Chem_m[i,j] < 0):
            b3[j] += s_eq*abs(Chem_m[i,j])*gs_u  * E_inh

    b = b1+b2+b3
    

    # linear system Av = b
    v_th = linalg.solve(A, transpose(b)) * volt

    return v_th
    
    
    
