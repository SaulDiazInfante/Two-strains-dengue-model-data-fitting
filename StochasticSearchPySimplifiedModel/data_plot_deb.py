import numpy as np
import datetime
import matplotlib.pyplot as plt
from stochastic_search import StochasticSearch
import os
sim = StochasticSearch()
sim.__init__()
sim.number_of_samples = 10
sim.bound_error_FD = 240
sim.bound_error_FHD = 50
##
sim.parameters_sampling_per_week()
sim.ode_int_solution_per_week()
sim.solution_plot()
#
t = sim.t
Y_m1_h = sim.solution[:, 8]
z = sim.solution[:, 12]
t_data_DF = sim.frecuency_per_week_DF[3:, 0]
t_data_DHF = sim.frecuency_per_week_DHF[1:, 0]
offset = 10000
#
t_z = t[0: -1: offset]
t_z = np.round(t_z)
t_z = t_z.astype(int)
z_points = z[0:-1: offset]

#
delete_index_t_z = [0, 1, 6, 9]
t_z = np.delete(t_z, delete_index_t_z)
z_points = np.delete(z_points, delete_index_t_z)
#
t_Y_m1_h = t[0: -1: offset]
t_Y_m1_h = np.round(t_Y_m1_h)
t_Y_m1_h = t_Y_m1_h.astype(int)
delte_index_t_Y = [0, 1, 3, 4, 6, 7, 8, 9]
t_Y_m1_h = np.delete(t_Y_m1_h, delte_index_t_Y)
Y_m1_h_points = Y_m1_h[0: -1: offset]
Y_m1_h_points = np.delete(Y_m1_h_points, delte_index_t_Y)

#
frecuency_per_week_DF = sim.frecuency_per_week_DF[3:, 1]
frecuency_per_week_DHF = sim.frecuency_per_week_DHF[1:, 1]
fitting_error_DF = \
    np.linalg.norm(frecuency_per_week_DF - z_points, ord=np.inf)
sim.fitting_error_DF = fitting_error_DF
fitting_error_DHF = \
    np.linalg.norm(frecuency_per_week_DHF - Y_m1_h_points, ord=np.inf)
sim.fitting_error_DHF = fitting_error_DHF
#
f1, ax_array = plt.subplots(2, 2, sharex=True)
ax_array[0, 0].plot(t, z, 'b-')
ax_array[0, 0].set_title(r'Reported DF ')

ax_array[0, 1].plot(t_data_DF, frecuency_per_week_DF,
		    ls='--',
		    color='lightblue',
		    marker='o',
		    ms=8,
		    mfc='lightblue',
		    alpha=0.7)
ax_array[0, 1].plot(t, z,
		    ls=':',
		    color='darkblue',
		    alpha=0.3)
ax_array[0, 1].plot(t_z, z_points,
		    ls='none',
		    color='blue',
		    marker='*',
		    ms=8,
		    mfc='blue',
		    alpha=0.5)
ax_array[0, 1].text(27, 300,
		    'err=' + str(np.round(fitting_error_DF, 1)),
		    fontsize=10
		    )
ax_array[0, 1].set_ylim(0, 400)

ax_array[0, 1].set_title(r'DF Fitting ')
#
ax_array[1, 0].plot(t, Y_m1_h, 'r-')
ax_array[1, 0].set_title(r'DHF')
ax_array[1, 1].plot(t_data_DHF, frecuency_per_week_DHF,
		    ls='--',
		    color='orange',
		    marker='o',
		    ms=8,
		    mfc='orange',
		    alpha=0.5)
ax_array[1, 1].plot(t, Y_m1_h,
		    ls=':',
		    color='crimson',
		    alpha=0.5
		    )
ax_array[1, 1].plot(t_Y_m1_h, Y_m1_h_points,
		    ls='none',
		    color='crimson',
		    marker='*'
		    )
ax_array[1, 1].text(27, 50,
		    'err=' + str(np.round(fitting_error_DHF, 1)),
		    fontsize=10
		    )
ax_array[1, 1].set_ylim(0, 60)
ax_array[1, 1].set_title(r'DHF Fitting ')

for i in np.arange(2):
    ax_array[1, i].set(xlabel='week n')
for j in np.arange(2):
    ax_array[j, 0].set(ylabel='Individuals')

plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
plt.show()