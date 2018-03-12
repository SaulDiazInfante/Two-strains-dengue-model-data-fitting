import numpy as np
import datetime
from stochastic_search import StochasticSearch
import os
sim = StochasticSearch()
sim.__init__()
sim.number_of_samples = 100000
sim.bound_error_FD = 240
sim.bound_error_FHD = 50
print '%-5s%-12s%-12s%-12s%-12s%-12s' \
      % ('i', 'R_01', 'R_02', 'R_zero', 'error_DF', 'error_DHF')
for i in np.arange(sim.number_of_samples):
    sim.parameters_sampling_per_week()
    sim.ode_int_solution_per_week()
    sim.solution_plot()
    sim.fitting_plot()
    error_DF = sim.fitting_error_DF
    error_DHF = sim.fitting_error_DHF
    r01_per_week, r02_per_week, r0_per_week = sim.compute_r_zero()
    sim.update_conditions_search()
    if sim.stop_condition:
        sim.save_parameters()
        str_time = str(datetime.datetime.now())
        file_name = './plots/fitting_DF_DHF' + str_time + '.png'
        os.rename('./plots/fitting_DF_DHF.png', file_name)
        file_name = './plots/populations_grid.png' + str_time + '.png'
        os.rename('./plots/populations_grid.png', file_name)
        print '=)'
    print '%-5d%-12f%-12f%-12f%-12f%-12f' \
          % (i, r01_per_week, r02_per_week, r0_per_week, error_DF, error_DHF)
