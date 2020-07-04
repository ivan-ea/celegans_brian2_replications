# -*- coding: utf-8 -*-
'''
Created on Sun May 25 10:03:20 2020

@author: Ivan Estevez Albuja
'''

from brian2 import * 
from aux_functions import  normalize_trace_avg, save_plot_list,normalize_trace_avg_unitary
from aux_functions import norm_trace_avg_unitary_2points
from Izq_model import Izq_gap_jun, Izq_chem_syn, Izq_eqs_constIn
#from genetic_functions import random_population, test_individual, build_matrices
from genetic_functions import fitness, build_matrices, test_individual, f1,f2,f3,f3_mod
import sys

N = 11 # number of neurons
time_step = 0.0025*second
defaultclock.dt = time_step

duration = 26*second;   transient_time = 6*second;  eval_time = 20*second
# duration = 60*second;   transient_time = 0*second;  eval_time = 20*second

n_points = int(duration/defaultclock.dt)
transient_points = int(transient_time/defaultclock.dt)
eval_points = int(eval_time/defaultclock.dt)
duration_time_series = arange(n_points)*time_step
eval_start = n_points - eval_points
eval_time_series = duration_time_series[:eval_points]
last_sec_points = int(2/(time_step/second))
last_sec_start = n_points - last_sec_points
last_sec_ev_start = eval_points - last_sec_points
s1= last_sec_start
s2= last_sec_start + n_points
e1= n_points
e2= 2*n_points        

labels = ['ASa', 'ASp', 'DAa', 'DAp', 'VAa', 'VAp', 'VBa', 'VBp', 'VDa', 'VDp', 'DB']
# index       0,     1,     2,     3,     4,     5,     6,     7,     8,     9,   10  

fwd_mask = array([0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1])
bkw_mask = array([0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0])

def get_trace(indv):
    neurons = NeuronGroup(N, model = Izq_eqs_constIn, method = 'euler')
    neurons.v = 0
    
    Chem_mat_test, Gap_mat_test = build_matrices(indv)
    
    # Chemical Synapses
    chem_sources, chem_targets = Chem_mat_test.nonzero()
    C = Synapses(neurons, neurons, model = Izq_chem_syn)
    C.connect(i = chem_sources, j = chem_targets)
    C.w = 0 
        
    # Gap Junctions
    gap_sources, gap_targets = Gap_mat_test.nonzero()
    G = Synapses(neurons, neurons, model = Izq_gap_jun)
    G.connect(i = gap_sources, j = gap_targets)
    G.g = 0
        
    trace = StateMonitor(neurons, 'v', record=True)
          
    Chem_mat, Gap_mat = build_matrices(indv)
    for i, j in zip(chem_sources, chem_targets):
        C.w[i,j] = Chem_mat[i,j]
    for i, j in zip(gap_sources, gap_targets):
        G.g[i,j] = Gap_mat[i,j]
    
    for i in range(0,10,2):
        neurons.bias[i] = neurons.bias[i+1] = indv[20 + int(i/2)]
        neurons.tau[i] = neurons.tau[i+1] = indv[26 + int(i/2)] * second
    neurons.bias[10] = indv[25]
    neurons.tau[10] = indv[31] * second  
        
    # forwards simulation
    neurons.active_ = fwd_mask
    neurons.I_ext = indv[32] * fwd_mask
    
    run(duration)
    
    # backwards simulation
    neurons.active_ = bkw_mask
    neurons.I_ext = indv[33] * bkw_mask
    
    run(duration)  
    
    # forwards simulation # only for final testing fw(10s) bkw(20s) fw(10s), no transient
    # neurons.active_ = fwd_mask
    # neurons.I_ext = indv[32] * fwd_mask
    # run(duration)  
    
    return trace.v


def fitfunc(indv):     
    trace = get_trace(indv)
    
    norm_trace = norm_trace_avg_unitary_2points(trace,s1,e1,s2,e2)
    

    fwd_trace=(norm_trace[:, eval_start:n_points])

    bkw_trace=(norm_trace[:, n_points+eval_start:2*n_points])
    
    # norm_f_trace = normalize_trace_avg(fwd_trace, eval_points-4*last_sec_points)
    # norm_b_trace = normalize_trace_avg(bkw_trace, eval_points-4*last_sec_points)
    
    # return fitness(norm_f_trace, norm_b_trace, time_step/second)
    return fitness(fwd_trace, bkw_trace, time_step/second)
  

def fitfunc_verbose(indv, save_traces = False):
    trace = get_trace(indv)
    
    norm_trace = norm_trace_avg_unitary_2points(trace,s1,e1,s2,e2)
    fwd_trace=(trace[:, eval_start:n_points])
    
    bkw_trace=(trace[:, n_points+eval_start:2*n_points])
    
    norm_f_trace = (norm_trace[:, eval_start:n_points])
    norm_b_trace = (norm_trace[:, n_points+eval_start:2*n_points])

    # plots 
    
    plot_labels = ['AS', 'DA', 'VA', 'DB', 'VB']
    plot_indices = [0,2,4,10,6]
    for j in plot_indices:
        plot(eval_time_series, fwd_trace[j], label = 'v_'+labels[j])
        xlabel('Time (s)')
        ylabel('V')
        title('Forward traces')
        legend(loc='upper right')
    show() 
    
    for j in plot_indices:
        plot(eval_time_series, norm_f_trace[j], label = 'v_'+labels[j])
        xlabel('Time (s)')
        ylabel('dV')
        title('Norm Forward traces')
        legend(loc='upper right')
    show() 
        
    for j in plot_indices:
        plot(eval_time_series, bkw_trace[j], label = 'v_'+labels[j])
        xlabel('Time (s)')
        ylabel('V')
        title('Backward traces')
        legend(loc='upper right')
    show() 
    
    for j in plot_indices:
        plot(eval_time_series, norm_b_trace[j], label = 'v_'+labels[j])
        xlabel('Time (s)')
        ylabel('dV')
        title('Norm Backward traces')
        legend(loc='upper right')
    show() 
    
   #only for final testing fw+bkw+fw 
    # for j in plot_indices:   
    #     plot(arange(3*n_points)*time_step ,trace[j], label = 'v_'+labels[j])    
    #     xlabel('Time (s)')
    #     ylabel('V')
    #     title('Full traces (forward + backward)')
    #     legend(loc='upper right')       
    # show() 
    
    # for j in plot_indices:   
    #     plot(arange(3*n_points-transient_points)*time_step ,norm_trace[j][:], label = 'v_'+labels[j])    
    #     xlabel('Time (s)')
    #     ylabel('V')
    #     title('Full traces norm (forward + backward)')
    #     legend(loc='upper right')       
    # show() 
    
    print('fitness normalized')
    print('f1=', f1(norm_f_trace, norm_b_trace, time_step/second, eval_time/second)) 
    print('f2=', f2(norm_f_trace, norm_b_trace, time_step/second, eval_time/second))
    print('f3=', f3(norm_f_trace, norm_b_trace))
    print('f3_mod=', f3_mod(norm_f_trace, norm_b_trace))
    f = fitness(norm_f_trace, norm_b_trace, time_step/second)
    print('f=', f)
    
    print('fitness not normalized')
    print('f1=', f1(fwd_trace, bkw_trace, time_step/second, eval_time/second)) 
    print('f2=', f2(fwd_trace, bkw_trace, time_step/second, eval_time/second))
    print('f3=', f3(fwd_trace, bkw_trace))
    print('f3_mod=', f3_mod(fwd_trace, bkw_trace))
    f = fitness(fwd_trace, bkw_trace, time_step/second)
    print('f=', f)
    
    
    if (save_traces):
        def save (name_trace, trace):
            name = 'I_Gen_'+name_trace+'_'+str(time_step/ms)+'_ev'+str(int(eval_time/second))
            head = 't (s)\t'        
            series_list = [duration_time_series[0:eval_points]]
            
            for i,ix in enumerate(plot_indices):
                head = head + plot_labels[i]+'\t' 
                series_list += [trace[ix]]
            
            save_plot_list(name, head, series_list)
            print("Written data in: ", name)
            return
        
        names = ['norm_Ftrace','F_trace','norm_Btrace','B_trace']
        traces = [norm_f_trace, fwd_trace, norm_b_trace, bkw_trace]
        
        for i in range(0,4):
            save(names[i],traces[i])
            
        name = 'I_Gen_'+'Full_trace_'+str(time_step/ms)+'_dur'+str(int(3*duration/second))
        head = 't (s)\t'        
        series_list = [arange(3*n_points)*time_step]
        
        for i,ix in enumerate(plot_indices):
            head = head + plot_labels[i]+'\t' 
            series_list += [trace[ix]]
        
        save_plot_list(name, head, series_list)
        print("Written data in: ", name)    
        
        
        name = 'I_Gen_'+'Full_trace_norm_'+str(time_step/ms)+'_dur'+str(int(3*duration/second))
        head = 't (s)\t'        
        series_list = [arange(3*n_points)*time_step]
        
        for i,ix in enumerate(plot_indices):
            head = head + plot_labels[i]+'\t' 
            series_list += [norm_trace[ix]]
        
        save_plot_list(name, head, series_list)
        print("Written data in: ", name)  
        
    
    return f   

