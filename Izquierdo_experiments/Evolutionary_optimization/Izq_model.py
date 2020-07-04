# -*- coding: utf-8 -*-
"""
Created on Mon May 11 11:31:45 2020

@author: Ivan Estevez Albuja

# Neuron model for C. elegans by Izquierdo and Beer.
# Original publication:  https://doi.org/10.1371/journal.pcbi.1002890
"""

from brian2 import * 


# Parameters

# Equations 
Izq_eqs_noInput = '''
dv/dt = (I_ext*active - v + Igap + Isyn) / tau : 1
Igap : 1 # gap junction current
Isyn : 1 # chemical synapsis current
# per neuron parameters
tau : second (constant)
bias : 1 (constant)
active: 1 (constant)
'''


Izq_eqs = Izq_eqs_noInput + 'I_ext = I_recorded(t) : 1' + ' \n'

Izq_eqs_constIn = Izq_eqs_noInput + 'I_ext : 1 (constant)' + ' \n'


Izq_gap_jun = '''
g : 1 # gap junction 'conductance', always > 0
Igap_post = g * (v_pre - v_post) : 1 (summed)
'''

Izq_chem_syn = '''
w : 1 # synaptic weight: >0 exhitatory, < 0 inhibitory
Isyn_post = w/(1+exp(-(v_pre + bias_pre))): 1 (summed)
'''
