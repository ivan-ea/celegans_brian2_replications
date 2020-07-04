# -*- coding: utf-8 -*-
'''
@author: Ivan Estevez Albuja
'''

from Brian2_experiment import brian2sim, labels, N, idx_labels_start_with
from Brian2_experiment import motor_indices, fwd_indices, fwd_indices_rel, all_fwd_indices
from Brian2_experiment import loco_indices, loco_colors, loco_indices_rel, all_loco_indices
from aux_functions import normalize_trace_avg, save_plot_list
from aux_functions import matrix_shape, plot_ntype_color, plot_traces_color
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 20})
import sys

#times in seconds
t_sim = 12
t_eval = 10
t_step = 0.00005

input_neuron = 'PLM'

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

rec_vars=['v']

in_indices = idx_labels_start_with(input_neuron, labels)

rec_indices = range(N)
# rec_indices = in_indices + motor_indices

trace = brian2sim(t_sim, t_step, I_max, rec_vars,
                  in_indices, rec_indices, abl_indices=None)


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

norm_mtraces = norm_v_traces[motor_indices,:]    
        
P, Sigma, Q_t = np.linalg.svd(norm_mtraces, full_matrices = False)

plt.plot(range(1,6), Sigma[0:5] * Sigma[0:5]/sum(Sigma*Sigma),'o')
plt.xlabel('n')
plt.ylabel('sig^2/sum(sig^2)')
plt.show()
print('# n \t sigma\t sig^2/sum(sig^2)')
for i in range(5):
    print(i+1,Sigma[i],Sigma[i]*Sigma[i]/sum(Sigma*Sigma))
# plt.plot(range(1,6), Sigma[0:5],'o')
# plt.show()

a1 = Sigma[0]*Q_t[0,:]
a2 = Sigma[1]*Q_t[1,:]

v1 = Sigma[0] * np.outer(P[:,0], Q_t[0,:])
v2 = Sigma[1] * np.outer(P[:,1], Q_t[1,:])

# comprobado que
#col = P[:,0]*a1[j] + P[:,1]*a2[j] - (v1+v2)[:,j]

plt.plot(eval_time_series, a1, label = 'a1')
plt.plot(eval_time_series, a2, label = 'a2')
plt.xlabel('Time (s)')
plt.ylabel('a')
plt.legend(loc='upper right')
plt.show()  

plt.plot(a1, a2)
plt.xlabel('a1')
plt.ylabel('a2')
plt.show()

sys.exit()

name = 'K_SVD_h'+str(t_step*1000)+'_I'+str(I_max)+'_ev'+str(int(t_eval))
head = 't (s)\tI (pA)\t a1 \t a2' 
series_list = [duration_time_series[0:eval_points], I_shape[0:eval_points], a1, a2]
save_plot_list(name, head, series_list)
print("Written data in: ", name)

# see decomposition of all fwd indices
all_fwd_indices_rel = [motor_indices.index(i) for i in all_fwd_indices]


fig = plt.figure(figsize=[10,12])
plt.imshow(norm_v_traces[all_loco_indices,:], aspect = 'auto')
plt.show()

matrix_shape(loco_indices, norm_v_traces)

plot_traces_color(fwd_indices_rel, v1, eval_time_series, loco_colors, 'v1 fwd neurons')
plot_traces_color(fwd_indices_rel, v2, eval_time_series, loco_colors, 'v2 fwd neurons')
plot_traces_color(fwd_indices_rel, v1+v2, eval_time_series, loco_colors, 'v1+v2 fwd neurons')
plot_traces_color(fwd_indices, norm_v_traces, eval_time_series, loco_colors,'norm fwd traces')

plot_traces_color(loco_indices_rel, v1, eval_time_series, loco_colors, 'v1 locomotion neurons')
plot_traces_color(loco_indices_rel, v2, eval_time_series, loco_colors, 'v2 locomotion neurons')
plot_traces_color(loco_indices_rel, v1+v2, eval_time_series, loco_colors, 'v1+v2 locomotion neurons')
plot_traces_color(loco_indices, norm_v_traces, eval_time_series, loco_colors,'norm loco traces')

sys.exit()
#columnas por ambos metodos son iguales!
def check_columns(eps):
    posibles = 0
    diffs = 0
    for j in range(0,200000,10):   
        col = P[:,0]*a1[j] + P[:,1]*a2[j] - (v1+v2)[:,j]  
        for val in col:
            if abs(val) > eps:
                posibles +=1
                diffs += abs(val)
                #print(val)
        if j%100 == 0: print('j',j,'posibles ',posibles, 'difference ', diffs)
    print('posibles ', posibles, 'difference ', diffs)
    return

