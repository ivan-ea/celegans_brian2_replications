# -*- coding: utf-8 -*-
'''
@author: Ivan Estevez Albuja
'''

from Brian2_experiment import brian2sim, labels, N, idx_labels_start_with
from Brian2_experiment import motor_indices, sensory_indices, inter_indices, fwd_indices, all_fwd_indices
from Brian2_experiment import loco_indices, loco_colors, loco_indices_rel, all_loco_indices
from aux_functions import normalize_trace_avg, save_plot_list
from aux_functions import matrix_shape, plot_ntype_color, plot_traces_color
import numpy as np
import matplotlib.pyplot as plt
import sys

#times in seconds
t_sim = 12
t_eval = 10
t_step = 0.00005

n_points = int(t_sim/t_step)
eval_points = int(t_eval/t_step)
eval_start = n_points - eval_points
last_sec_points = int(1/t_step)
last_sec_start = n_points - last_sec_points
last_sec_ev_start = eval_points - last_sec_points

duration_time_series = np.arange(n_points)*t_step
eval_time_series = duration_time_series[:eval_points]

# External input waveform
I_max = 2000 # CHANGE THIS to 1000 and 1500

I_shape = I_max * np.ones(n_points)
#no oscillations at I = 1000 # oscillations start at I = 1500
#good at I=2000, weird at I=5000, stop at I=10000

rec_vars=['v']
in_indices = idx_labels_start_with('PLM', labels)

trace = brian2sim(t_sim, t_step, I_max, rec_vars,
                         in_indices, rec_indices = True, abl_indices=None)


# Plot Input without recording I_ext variable
plt.plot(duration_time_series, I_shape)
plt.xlabel('Time (s)')
plt.ylabel('I (pA)')
plt.show() 

v_traces = np.zeros([N, eval_points])
for i in range(N):
    v_traces[i] = trace[i].v[eval_start:]

norm_v_traces = normalize_trace_avg(v_traces, eval_points-last_sec_points)   
print('end normalization')

# average of traces by neuron type
types_indices = {'m': motor_indices, 's': sensory_indices, 'i': inter_indices}
types_averages = {'m': np.zeros(eval_points), 's': np.zeros(eval_points), 'i': np.zeros(eval_points)}   
for n_type in list(types_averages):
    for i in types_indices[n_type]:
        types_averages[n_type] = types_averages[n_type] + norm_v_traces[i]
    types_averages[n_type] = types_averages[n_type]/len(types_indices[n_type])  
print('end averages')

#plot average

for i in [0,1,2]:
    plt.plot(eval_time_series, list(types_averages.values())[i]*1000, 
         label = list(types_averages)[i])
    plt.title('averages for I = '+str(I_max))
    plt.xlabel('Time (s)')
    plt.ylabel('Delta V (mV)')
    plt.legend(loc='upper right')
plt.show()

matrix_shape(fwd_indices, norm_v_traces)
plot_traces_color(fwd_indices, norm_v_traces, eval_time_series, loco_colors,'norm fwd traces')

sys.exit()

name = 'Kunert_bif/K_bif_avgs_h'+str(t_step*1000)+'_I'+str(I_max)+'_eval'+str(t_eval)
head = 't (s)\t d_avg_m (mV)\td_avg_s (mV)\td_avg_i (mV)'
s = [duration_time_series[0:eval_points]]
for t in ['m','s','i']:
    s.append(types_averages[t][:]*1000)
save_plot_list(name, head, s)
print("Written data in: ", name)    


name = 'Kunert_bif/K_bif_MotorMatrix'+str(t_step*1000)+'_I'+str(I_max)+'_ev'+str(t_eval)
head = ''
series_list = []
for i in all_fwd_indices:
    head = head + '\tdV_'+labels[i]+' (mv)'
    series_list.append(norm_v_traces[i]*1000)
save_plot_list(name, head, series_list)
print("Written data in: ", name)
