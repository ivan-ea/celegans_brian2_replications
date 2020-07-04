# -*- coding: utf-8 -*-
"""
@author: Ivan Estevez Albuja
"""

from brian2 import * 

#lineal triangle from y = -30 to y = + 30
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

def save_plot(filename,header,*series):
    file = open(filename + ".txt",'w')
    file.write("# "+header + "\n")
    n = len(series[0])
    for i in range(0,n):
        for serie in series:
            file.write("%0.3f"%(serie[i]) + " ")
        file.write("\n")
    return 0

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


