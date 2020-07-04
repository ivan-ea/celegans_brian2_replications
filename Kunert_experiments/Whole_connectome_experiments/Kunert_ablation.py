# -*- coding: utf-8 -*-
'''
@author: Ivan Estevez Albuja
'''

from Brian2_experiment import brian2sim, labels, N, idx_labels_start_with
from Brian2_experiment import motor_indices, fwd_indices, fwd_indices_rel, all_fwd_indices
from Brian2_experiment import loco_indices, loco_colors, loco_indices_rel, all_loco_indices
from aux_functions import normalize_trace_avg, save_plot_list, separated_plots_color
from aux_functions import matrix_shape, separated_plots_color, plot_traces_color
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 20})
import sys

#times in seconds
t_sim = 5
t_eval = 2
t_step = 0.00005

abl_runs = {
    'None': {'ablate': False},
    'AVA': {'ablate': True}, 
    'AVB': {'ablate': True}, 
    'AIZR': {'ablate': True}
    }


n_points = int(t_sim/t_step)  
eval_points = int(t_eval/t_step)
eval_start = n_points - eval_points
last_sec_points = int(1/t_step)
last_sec_start = n_points - last_sec_points
last_sec_ev_start = eval_points - last_sec_points

duration_time_series = np.arange(n_points)*t_step
eval_time_series = duration_time_series[:eval_points]

# External input waveform
I_max = 2000
I_shape = I_max * np.ones(n_points)
# Plot Input without recording I_ext variable
# plt.plot(duration_time_series, I_shape)
# plt.xlabel('Time (s)')
# plt.ylabel('I (pA)')
# plt.show()


rec_vars=['v']

in_indices = in_indices = idx_labels_start_with('PLM', labels)
    
rec_indices = range(N)

for ab in abl_runs:
    ablation_type = ab
    ablate = abl_runs[ab]['ablate']
    
    if ablate: 
        abl_indices = idx_labels_start_with(ablation_type, labels)
    else:
        abl_indices = None
    
    
    trace = brian2sim(t_sim, t_step, I_max, rec_vars,
                      in_indices, rec_indices, abl_indices)

    #Results after ablation 
    v_traces = np.zeros([N, eval_points])
    for i in range(N):
        v_traces[i] = trace[i].v[eval_start:]
    
    norm_v_traces = normalize_trace_avg(v_traces, eval_points-last_sec_points)   
    print('end normalization','ablation of',ab)
     
    norm_mtraces = norm_v_traces[motor_indices,:]      
            
    P, Sigma, Q_t = np.linalg.svd(norm_mtraces, full_matrices = False)
    
    abl_runs[ab]['sigma'] = (Sigma[0:5]*Sigma[0:5])/sum(Sigma*Sigma)

    print('# n \t sigma\t sig^2/sum(sig^2)','ablation of',ab)
    for i in range(5):
        print(i+1,Sigma[i],Sigma[i]*Sigma[i]/sum(Sigma*Sigma))

# plt.plot(range(1,6), Sigma[0:5],'o')
# plt.show()

    abl_runs[ab]['a1'] = Sigma[0] * Q_t[0,:]
    abl_runs[ab]['a2'] = Sigma[1] * Q_t[1,:]

    abl_runs[ab]['v1'] = Sigma[0] * np.outer(P[:,0], Q_t[0,:])
    abl_runs[ab]['v2'] = Sigma[1] * np.outer(P[:,1], Q_t[1,:])
    
    separated_plots_color(loco_indices, norm_v_traces, eval_time_series, loco_colors,'norm loco traces with '+ablation_type+' ablated')
    matrix_shape(loco_indices, norm_v_traces)

fig = plt.figure()
ax = fig.add_subplot(111)
for ab in abl_runs:
    ax.plot(abl_runs[ab]['a1'], abl_runs[ab]['a2'], label = ab)
plt.xlabel('a1')
plt.ylabel('a2')
plt.legend(loc='upper right')
ax.set_aspect('equal')
plt.show()

data = [abl_runs[ab]['sigma'][0:2] for ab in abl_runs]

X = np.arange(2)
fig = plt.figure()
ax = fig.add_axes([0,0,1,1])
ax.bar(X + 0.00, data[0], color = 'black', width = 0.25, label = "none")
ax.bar(X + 0.25, data[1], color = 'b', width = 0.25, label = "AVA")
ax.bar(X + 0.50, data[2], color = 'g', width = 0.25, label = "AVB")
ax.bar(X + 0.75, data[3], color = 'r', width = 0.25, label = "AIZR")
ax.legend(loc='upper right')


sys.exit()
# plt.plot(eval_time_series, a1, label = 'a1')
# plt.plot(eval_time_series, a2, label = 'a2')
# plt.xlabel('Time (s)')
# plt.ylabel('a')
# plt.legend(loc='upper right')
# plt.title(ablation_type+' ablations')
# plt.show()  

# plt.plot(a1, a2)
# plt.xlabel('a1')
# plt.ylabel('a2')
# plt.show()

# see decomposition of all fwd indices


# plot_traces_color(loco_indices, norm_v_traces, eval_time_series, loco_colors,'norm loco traces with '+ablation_type+' ablated')


#
sys.exit()

name = 'K_ABL_a1_a2_h'+str(t_step*1000)+'_I'+str(I_max)+'_ev'+str(t_eval)
head = 't (s)\tI (pA)'

series_list = [duration_time_series[0:eval_points], I_shape[0:eval_points]]
for ab in abl_runs:
    head = head + '\t a1' +ab + '\t a2' +ab 
    series_list += [abl_runs[ab]['a1'],abl_runs[ab]['a2']]

save_plot_list(name, head, series_list)
print("Written data in: ", name)









