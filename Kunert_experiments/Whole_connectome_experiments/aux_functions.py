# -*- coding: utf-8 -*-
"""
@author: Ivan Estevez Albuja
"""

from brian2 import * 

#lineal triangle from y = -Imax to y = + Imax
def triang_shape(t_rec_arr,delay,duration,Imax): 
    I_arr = zeros(size(t_rec_arr))
    nu_dur = duration - 2*delay*ms
    for i,t_rec in enumerate(t_rec_arr):
        if (t_rec < nu_dur/2):
            I_arr[i]=(2*Imax)*t_rec/(nu_dur/2)-Imax
        else:
            I_arr[i]=(-2*Imax)*t_rec/(nu_dur/2)+3*Imax       
    return I_arr

def trimmed_triangle(delay, duration,Imax):
    delayed_t = arange(int((duration-(2*delay)*ms)/defaultclock.dt))*defaultclock.dt
    shape = triang_shape(delayed_t,delay, duration,Imax)
    trim_t = arange(int((delay*ms)/defaultclock.dt))*defaultclock.dt
    I_arr = zeros(size(trim_t))
    I_arr[:] = -Imax
    return concatenate([I_arr,shape,I_arr])
    
def trimmed_rectangle(delay, Imax, duration):
    trim_t = arange(int((delay*ms)/defaultclock.dt))*defaultclock.dt
    I_arr = zeros(size(trim_t))

    shape = zeros(int((duration-(2*delay)*ms)/defaultclock.dt))
    shape[:] = Imax
    return concatenate([I_arr,shape,I_arr])

def trimmed_rectangle_tail(delay, Imax, duration):
    trim_t = arange(int((delay*ms)/defaultclock.dt))*defaultclock.dt
    I_arr = zeros(size(trim_t))

    shape = zeros(int((duration-(4*delay)*ms)/defaultclock.dt))
    shape[:] = Imax
    return concatenate([I_arr,shape,I_arr,I_arr,I_arr])

def square_wave(Imax,n_points):
    div_size = int(n_points/4)
    zeros_ = zeros(div_size)
    ones_ = ones(div_size)*Imax
    return concatenate([ones_,zeros_,ones_,zeros_])

def double_sqwave(Imax, n_points):
    sqwave = square_wave(Imax,int(n_points/2))
    return (concatenate([sqwave,sqwave]))

def n_sqwaves(n,Imax, n_points):
    sqwave = square_wave(Imax,int(n_points/n))
    res = concatenate([sqwave,sqwave])
    for i in range(2,n):
        res = concatenate([res,sqwave])
    return res

def constant_intput(Imax,duration):
    return ones(int(duration/defaultclock.dt))*Imax

'''time series to save are in a list'''
def save_plot_list(filename,header,series):
    file = open(filename + ".txt",'w')
    file.write("#" + header + "\n")
    n = len(series[0])
    for i in range(0,n):
        for serie in series:
            file.write("%0.7f"%(serie[i]) + " ")
        file.write("\n")
    return 0

'''time series to save as comma separated args'''
def save_plot(filename,header,*series):
    save_plot_list(filename,header,series)
    return 0


def diff_V(v_ts, E_ref):
    n = len(v_ts)
    return v_ts - ones(n)*E_ref

def normalize_trace_min(tr, start_point):
    norm_traces = zeros(shape(tr))
    for i in range(shape(tr)[0]):
        norm_traces[i] = diff_V(tr[i], min(tr[i][start_point:]))
    return norm_traces 

def normalize_trace_avg(tr, start_point):
    norm_traces = zeros(shape(tr))
    for i in range(shape(tr)[0]):
        mi = min(tr[i][start_point:])
        ma = max(tr[i][start_point:])
        norm_traces[i] = diff_V(tr[i], 0.5*(ma+mi))        
    return norm_traces 


def visualise_connectivity(S):
    Ns = len(S.source)
    Nt = len(S.target)
    figure(figsize=(10, 4))
    subplot(121)
    plot(zeros(Ns), arange(Ns), 'ok', ms=10)
    plot(ones(Nt), arange(Nt), 'ok', ms=10)
    for i, j in zip(S.i, S.j):
        plot([0, 1], [i, j], '-k')
    xticks([0, 1], ['Source', 'Target'])
    ylabel('Neuron index')
    xlim(-0.1, 1.1)
    ylim(-1, max(Ns, Nt))
    subplot(122)
    plot(S.i, S.j, 'ok')
    xlim(-1, Ns)
    ylim(-1, Nt)
    xlabel('Source neuron index')
    ylabel('Target neuron index')
    show()
    return

def matrix_shape(dic_indices,matrix):   
    n = len(dic_indices)
    fig, ax = plt.subplots(n, 1, sharex='col', sharey='row', figsize = [7,12])
    for i,label in enumerate(dic_indices):
        ax[i].imshow(matrix[dic_indices[label],:], aspect = 'auto')
    plt.show()
    return

def separated_plots_color(dic_indices, vtrace, times, colors, title):
    n = len(dic_indices)
    fig, ax = plt.subplots(n, 1, sharex='col', sharey='row', figsize = [7,12])
    for i,label in enumerate(dic_indices):
        for j in dic_indices[label]:
            ax[i].plot(times, vtrace[j]*1000, color = colors[label])
    #plt.ylabel('dV (mV)')
    # plt.xlabel('t (s)')
    plt.show()
    return

def plot_ntype_color(dic_indices, label, vtrace, times, colors, is_show):
    for i in dic_indices[label]:
        plt.plot(times, vtrace[i]*1000, color = colors[label])
    if (is_show):
        plt.ylabel('dV (mV)')
        plt.show()
    return


def plot_traces_color(dic_indices, vtrace, times, colors, title):
    for label in dic_indices:
        plot_ntype_color(dic_indices, label, vtrace, times, colors, False)
    plt.ylabel('dV (mV)')
    plt.title(title)
    plt.show()   
    return


