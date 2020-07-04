# -*- coding: utf-8 -*-
"""
Created on Thu Jun 11 21:37:14 2020

@author: Ivan Estevez Albuja
"""

import numpy as np
from geneticalgorithm import geneticalgorithm as ga 
#info # https://pypi.org/project/geneticalgorithm/

from Izq_CPG_fitness_func import fitfunc, fitfunc_verbose
from genetic_functions import tau_idx, g_idx

def f(x):
    return -1*fitfunc(x)

algorithm_param = {'max_num_iteration': 10,\
                   'population_size': 10,\
                   'mutation_probability': 0.3,\
                   'elit_ratio': 0.02,\
                   'crossover_probability': 0.5,\
                   'parents_portion': 0.25,\
                   'crossover_type':'uniform',\
                   'max_iteration_without_improv': 15}
    

varbound=np.array([[-20.0, 20.0]]*34)
for i in tau_idx:
    varbound[i] = [0.05, 2.0]
for i in g_idx:
    varbound[i] = [0.0, 2.5]



model=ga(function=f,  function_timeout=100, dimension=34,\
            variable_type='real',\
            variable_boundaries=varbound,\
            algorithm_parameters=algorithm_param)

model.run()

solution = model.best_variable

sol_fitness = model.best_function

fitfunc_verbose(solution)

print('#----------------------------------------------------------------------------')

best_so_far = np.array([-17.26358619,  -7.97575289,  -7.25144085,   3.57326048,
        13.38905508,   0.95439016, -19.76550766,  -9.21040112,
        10.17130201,  -6.57606904,   0.02832742,   1.14672932,
         1.0370267 ,   2.00993874,   0.16106674,  -5.52194449,
        13.84803218,   7.09557621,  -0.30618229,   5.20471193,
        -6.15974197,   7.18217651,  -1.41068216, -13.62878656,
         8.27304127, -13.05380308,   0.71525749,   0.71186105,
         0.53407642,   1.15593008,   1.90090529,   0.18609726,
         5.96782138,   0.92583596])


fitfunc_verbose(best_so_far)
