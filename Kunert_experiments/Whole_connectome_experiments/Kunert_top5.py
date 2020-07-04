# -*- coding: utf-8 -*-
'''
@author: Ivan Estevez Albuja
'''

from Brian2_experiment import brian2sim, labels, N, idx_labels_start_with
from aux_functions import n_sqwaves, save_plot_list, normalize_trace_min
import numpy as np
import matplotlib.pyplot as plt
import sys

#times in seconds
t_step = 0.00005

# isolated = True
# t_sim = 1
# t_eval = 1

isolated = False
t_sim = 5
t_eval = 1


n_points = int(t_sim/t_step)
eval_points = int(t_eval/t_step)
eval_start = n_points - eval_points
last_sec_points = int(1/t_step)
last_sec_start = n_points - last_sec_points
last_sec_ev_start = eval_points - last_sec_points

duration_time_series = np.arange(n_points)*t_step
eval_time_series = duration_time_series[eval_start:]

# External input waveform
I_max = 1
#I_shape = I_max * np.ones(n_points) # I constant
I_shape = n_sqwaves(t_sim,I_max,n_points) #

rec_vars = ['v', 's']
in_indices = idx_labels_start_with('ASHR', labels)
rec_indices = range(N) if (not isolated) else in_indices
abl_indices = None if (not isolated) else range(N)


trace = brian2sim(t_sim, t_step, I_shape, rec_vars, 
                  in_indices, rec_indices, abl_indices)


# Plot Input without recording I_ext variable
plt.plot(duration_time_series[0:eval_points], I_shape[0:eval_points])
plt.xlabel('Time (s)')
plt.ylabel('I (pA)')
plt.show()

vtraces = np.zeros([N, eval_points])
for i, ix in enumerate(rec_indices):
    vtraces[i] = trace[ix].v[eval_start:]  

norm_vtraces = normalize_trace_min(vtraces, eval_points-last_sec_points)   
print('end normalization')

if isolated:
    plt.plot(duration_time_series[:eval_points], norm_vtraces[0]*1000,label = 'V_'+labels[80])
    plt.xlabel('Time (s)')
    plt.ylabel('V (mV)')
    plt.show()
    sys.exit()
    name = 'Kunert_top5/K_top5_ISO_h'+str(t_step*1000)+'_I'+str(I_max)+'_eval'+str(t_eval)
    head = 't (s)\tI (pA)\tdV_ASHR (mv)\ts_ASHR' 
    series_list = [duration_time_series[0:eval_points], I_shape[0:eval_points], norm_vtraces[0]*1000]
    series_list.append(trace[80].s[eval_start:])
    save_plot_list(name, head, series_list)
    print("Written data in: ", name)
    sys.exit()

def calc_top5_vtraces(norm_vtraces):
    maxes = [] 
    for v_ in norm_vtraces:
        maxes.append(max(v_))        
    top_indices = (np.argsort(maxes))[::-1][0:6]  
    return top_indices
# top_indices = calc_top5_vtraces(norm_vtraces)
    
def calc_min5_vtraces(norm_vtraces, eval_points):    
    safety_pts = 2250 # a number of time steps for safety interval
    start_valley = int(3*eval_points/4) + safety_pts # 
    end_valley = int(eval_points) - safety_pts
    
    maxes_valley = []
    for v_ in norm_vtraces:
        maxes_valley.append(max(v_[start_valley:end_valley]))
    
    top_indices_valley = (np.argsort(maxes_valley))[::-1][0:5]   
    return top_indices_valley
# top_indices_valley = calc_min5_vtraces(norm_vtraces, eval_points)

top_indices = [80, 131, 157, 57, 76, 153]  # for Imax = 1
top_indices_valley = [134, 266, 121, 113, 69] # before, last was = 1
plot_indices = np.concatenate([top_indices[1:], top_indices_valley])

#V of ASHR connected
plt.plot(duration_time_series[:eval_points], norm_vtraces[80]*1000,label = 'V_'+labels[80])
plt.xlabel('Time (s)')
plt.ylabel('V (mV)')
plt.show()

#Plot V 
for i in plot_indices:    
    plt.plot(duration_time_series[:eval_points], norm_vtraces[i]*1000, 
         label = 'V_'+labels[i],)
    plt.xlabel('Time (s)')
    plt.ylabel('Delta V (mV)')
    plt.legend(loc='upper right')
plt.show()

k=0
colors = ['red','red','red','red','red','lime','lime','lime','lime','lime']
#Plot V 2 colors
for i in plot_indices:    
    plt.plot(duration_time_series[:eval_points], norm_vtraces[i]*1000, 
         label = 'V_'+labels[i], color = colors[k])
    plt.xlabel('Time (s)')
    plt.ylabel('Delta V (mV)')
    plt.legend(loc='upper right')
    k +=1
plt.show()

sys.exit()

#plot s
for i in plot_indices:    
    plt.plot(duration_time_series[:eval_points], trace[i].s[eval_start:],label='s_'+labels[i])
    plt.xlabel('Time (s)')
    plt.ylabel('s')
    plt.legend(loc='upper right')
plt.show()

#not normalized traces
k=0
for i in plot_indices:    
    plt.plot(duration_time_series[:eval_points], vtraces[i],
         label = 'V_'+labels[i], color = colors[k])
    plt.xlabel('Time (s)')
    plt.ylabel('Delta V (mV)')
    plt.legend(loc='upper right')
    k +=1
plt.show()

#original traces

for i in plot_indices:    
    plt.plot(duration_time_series[:], trace[i].v[:],label = 'V_'+labels[i])
    plt.xlabel('Time (s)')
    plt.ylabel('V (mV)')
    plt.legend(loc='upper right')
plt.show()

sys.exit()


name = 'Kunert_top5/K_top5_h'+str(t_step*1000)+'_I'+str(I_max)+'_eval'+str(t_eval)
head = 't (s)\tI (pA)\tdV_ASHR (mv)\ts_ASHR' 
series_list = [duration_time_series[0:eval_points], I_shape[0:eval_points], norm_vtraces[80]*1000]
series_list.append(trace[80].s[eval_start:])

for i in plot_indices:
    head = head + '\tdV_'+labels[i]+' (mv)'
    series_list.append(norm_vtraces[i]*1000)
    head = head + '\ts_'+labels[i]
    series_list.append(trace[i].s[eval_start:])

save_plot_list(name, head, series_list)
print("Written data in: ", name)

