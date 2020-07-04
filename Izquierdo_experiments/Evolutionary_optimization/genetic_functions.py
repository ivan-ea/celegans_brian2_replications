# -*- coding: utf-8 -*-
"""
Created on Mon May 25 12:11:33 2020

@author: Ivan Estevez Albuja
"""

''' individual = numpy array of 34 elements
    0-9 w_ij: chemical synapses weight
    10-13 : g_ij gap junction conductance
    14-31: neuron parameters w_ii (14-19), bias_i (20-25), tau_i (26-31)
    32, 33: I_ext from AVA and AVB
'''
labels = ['ASa', 'ASp', 'DAa', 'DAp', 'VAa', 'VAp', 'VBa', 'VBp', 'VDa', 'VDp', 'DB']
#             0,     1,     2,     3,     4,     5,     6,     7,     8,     9,   10 
import numpy as np

from numpy.random import rand
from scipy import integrate
import math

n_genes = 34
w_idx = np.concatenate([np.arange(0,10),np.arange(14,20)])
bias_idx = np.arange(20,26)
tau_idx = np.arange(26,32)
g_idx = np.arange(10,14)
i_idx = np.arange(32,34)

Amp = 0.3 #optimal oscillation amplitude
 
test_individual = np.arange(1,35)

def build_matrices(ind):
    N = 11
    
    Chem_mat = np.zeros([N,N])    
    Gap_mat = np.zeros([N,N])  
    diagonal = np.eye(N)
    
    # Chemical Synapses    
    # from AS
    Chem_mat[0,2] = Chem_mat[1,3] = ind[0]
    Chem_mat[0,8] = Chem_mat[1,9] = ind[1]   
    # from DA
    Chem_mat[2,8] = Chem_mat[3,9] = ind[2]
    Chem_mat[2,9] = Chem_mat[3,8] = ind[3]
    Chem_mat[2,10] = Chem_mat[3,10] = ind[4]    
    # from VA
    Chem_mat[4,8] = Chem_mat[5,9] = ind[5]    
    # from VD
    Chem_mat[8,4] = Chem_mat[9,5] = ind[6]
    Chem_mat[8,6] = Chem_mat[9,7] = ind[7]    
    # from DB
    Chem_mat[10,0] = Chem_mat[10,1] = ind[8]
    Chem_mat[10,8] = Chem_mat[10,9] = ind[9]

    for i in range(N-1):
        diagonal[i] *= ind[int(i/2)+14]
    diagonal[-1] *= ind[19] 
    
    Chem_mat += diagonal
    
    #Gap Junctions
    # DA -- VA
    Gap_mat[2,4] = Gap_mat[3,5] = ind[10]
    Gap_mat[4,2] = Gap_mat[5,3] = ind[10]
    # VD -- VA
    Gap_mat[8,4] = Gap_mat[9,5] = ind[11]
    Gap_mat[4,8] = Gap_mat[5,9] = ind[11]
    # VBa -- VBp
    Gap_mat[6,7] = Gap_mat[7,6] = ind[12]
    # VDa -- VDp
    Gap_mat[8,9] = Gap_mat[9,8] = ind[13] 
    
    return Chem_mat, Gap_mat

def scale_rand(min_, max_):
    r = rand()
    return min_ + (r * (max_ - min_))

def random_individual(n_genes):
    ind = np.zeros(n_genes)
    for i in range(n_genes):
        ind[i] = scale_rand(-20,20)
    for i in tau_idx:
        ind[i] = scale_rand(0.05, 2)
    for i in g_idx:
        ind[i] = scale_rand(0,2.5)
               
    return ind
    
def random_population(n):
    pop = np.zeros([n,n_genes])
    for i in range(n):
        pop[i] = random_individual(n_genes)
    return pop

    

def f1(f_tr, b_tr, dt, T): # Oscillation Criterion
    f_cells = ['DB', 'VBa','VBp','ASa']
    b_cells = ['DAa', 'DAp', 'VAa', 'VAp','ASp']
    f_indx = list(map(lambda x: labels.index(x), f_cells))    
    b_indx = list(map(lambda x: labels.index(x), b_cells))

    s = {}
    
    for i, cell in zip(f_indx, f_cells):
        I = integrate.trapz(y = abs(np.diff(f_tr[i])/dt), dx = dt)
        s[cell] = 2*I/(Amp*T)
        
    for i, cell in zip(b_indx, b_cells):
        I = integrate.trapz(y = abs(np.diff(b_tr[i])/dt), dx = dt)
        s[cell] = 2*I/(Amp*T)   
    
    #print('f1, s',s) #debug
    
    prod = 1
    for v in s.values():
        prod *= v
        
    #print('f1, prod', prod) #debug
    return prod

def f2(f_tr, b_tr, dt, T): # Phase Criterion
    f_set = [('VBa', 'DB'), ('VBp', 'DB')]
    b_set = [('VAa', 'DAa'), ('VAp','DAp')]
    
    f_phase_lag = 1   
    for v,d in f_set:
        v_series = np.diff(f_tr[labels.index(v)])
        d_series = np.diff(f_tr[labels.index(d)])
        I = integrate.trapz(y = abs(np.sign(v_series) + np.sign(d_series)), dx = dt)
        f_phase_lag = f_phase_lag * (1-I/(2*T))
    
    b_phase_lag = 1
    for v,d in b_set:
        v_series = np.diff(b_tr[labels.index(v)])
        d_series = np.diff(b_tr[labels.index(d)])
        I = integrate.trapz(y = abs(np.sign(v_series) + np.sign(d_series)), dx = dt)
        b_phase_lag = b_phase_lag * (1-I/(2*T))  
        
    #print('phase_lag:',f_phase_lag,b_phase_lag)  #debug  
    return f_phase_lag * b_phase_lag

def f3(f_tr, b_tr): # Dominance Criterion
    
    def smooth_fn(x, x_0):
        r = x/x_0
        f = 0.1 + 0.9 * r * math.exp(1-r)
        if f < 0: f = 0
        #print('r',r,' f_smooth', 0.1 + 0.9 * r * math.exp(1-r))#debug
        return f
    
    y = ['DB', 'VBa','VBp']
    x = ['DAa', 'DAp', 'VAa', 'VAp']
    
    dic = {}
    for key in y+x:
        dic[key]= {'m_y': 0, 'M_y': 0, 'M_x': 0}
        
    for l in y:
        dic[l]['m_y'] = min(f_tr[labels.index(l)])
        dic[l]['M_y'] = max(f_tr[labels.index(l)])
        dic[l]['M_x'] = max(b_tr[labels.index(l)])
        
    for l in x:
        dic[l]['m_y'] = min(b_tr[labels.index(l)])
        dic[l]['M_y'] = max(b_tr[labels.index(l)])
        dic[l]['M_x'] = max(f_tr[labels.index(l)])        
        
    prod_m_y = 1
    prod_A_y = 1
    prod_M_x = 1
    
    for l in dic:
        prod_m_y *= smooth_fn(dic[l]['m_y'], 1-Amp)
        prod_M_x *= smooth_fn(dic[l]['M_x'], Amp)
        prod_A_y *= smooth_fn(dic[l]['M_y']-dic[l]['m_y'], Amp)
        
    #print('f3:', '%0.4g'%prod_m_y, '%0.4g'%prod_A_y, '%0.4g'%prod_M_x) # debug   
    # min_ = min(prod_m_y, prod_A_y, prod_M_x, prod_m_y * prod_A_y * prod_M_x)
        
    return prod_m_y * prod_A_y * prod_M_x

def f3_mod(f_tr,b_tr):
    y = ['DB', 'VBa']#['DB', 'VBa','VBp']
    x = ['DAa','VAa']#['DAa', 'DAp', 'VAa', 'VAp']
    
    dic = {}
    for key in y+x:
        dic[key]= { 'M_y': 0, 'M_x': 0}
        
    for l in y:
        dic[l]['M_y'] = max(f_tr[labels.index(l)])
        dic[l]['M_x'] = max(b_tr[labels.index(l)])
        
    for l in x:
        dic[l]['M_y'] = max(b_tr[labels.index(l)])
        dic[l]['M_x'] = max(f_tr[labels.index(l)]) 
    
    score = 0.0

    for l in dic:
        if dic[l]['M_y']> dic[l]['M_x']:  score += 1
        
    for l in ['VBa','VAa','DAa']:
        if dic[l]['M_y'] > 0 and dic[l]['M_x'] < 0: score +=1
    
    return score
    

def fitness(f_trace, b_trace, dt):
    eval_points = np.shape(f_trace)[1]
    T = eval_points * dt
    # print('T=',T, 'dt=',dt) #debug
    f_1 = 1
    f_2 = 1
    f_3 = 1
    
    f_1 = f1(f_trace, b_trace,dt,T)
    if f_1 <= 1e-60: return 0
    # magnitude = math.floor(math.log10(f_1))
    # factor = 10** magnitude
    #print('magnitude',magnitude, 'factor',factor)    

    if f_1 < 1000000: return f_1*0.0000001

    f_3 = f3_mod(f_trace, b_trace)
    f_2 = f2(f_trace, b_trace,dt,T)

    f = 22*f_2 + 4*f_3 + 0.9*math.log10(f_1) #+ f_1/factor 
    # s = '%0.4g'
    # print('fitness: f1=',s%(f_1),', f2=',s%(10*f_2),', f3=',s%(f_3),', f=',s%f)
    return  f


