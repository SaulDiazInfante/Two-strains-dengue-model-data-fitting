from model_numerics_per_week import NumericOdeSolutionPerWeek
import numpy as np
import matplotlib.pyplot as plt


class StochasticSearch(NumericOdeSolutionPerWeek):
    def __init__(self):
        self.number_of_samples = 30
        self.frecuency_per_week_DF = \
            np.loadtxt('./data/frecuency_per_week_DF.dat', dtype='int',
                       delimiter=',')
        self.frecuency_per_week_DHF = \
            np.loadtxt('./data/frecuency_per_week_DHF.dat', dtype='int',
                       delimiter=',')
        self.bound_error_FD = 100
        self.bound_error_FHD = 30
        self.bound_initial_error_FHD = 20
        self.bound_initial_error_FD = 10
        # Numerics initial parameters
        self.Lambda_M = 7 * 20000.000000
        self.Lambda_S = 7 * 1.400000
        self.Lambda_m1 = 7 * 0.015000
        self.Lambda_m2 = 7 * 0.065000
        self.beta_M = 0.029149
        self.beta_H = 0.030855
        self.b = 7 * 2.095346
        self.mu_M = 7 * 0.076923
        self.mu_H = 7 * 0.000039
        self.alpha_c = 7 * 0.256219
        self.alpha_h = 7 * 0.256219
        self.sigma = 5.000000
        self.p = 0.050000
        self.theta = 0.450000
        self.M_s0 = 120000.000000
        self.M_10 = 20.000000
        self.M_20 = 30.000000
        self.I_s0 = 35600.000000
        self.I_10 = 1.000000
        self.I_20 = 20.000000
        self.S_m10 = 4400.000000
        self.S_m20 = 0.0
        self.Y_m1_c0 = 0.0
        self.Y_m1_h0 = 0.0
        self.Y_m2_c0 = 0.0
        self.Rec_0 = 0.0
        self.z0 = 1.050000
        self.t0 = 25  # 25.0
        self.T = 53
        self.grid_size = np.int(self.T - self.t0) * 10000
        self.h = np.float64(self.T) / np.float64(self.grid_size)
        self.r_01 = 0.0
        self.r_02 = 0.0
        self.r_zero = 0
        #
        self.t = np.linspace(self.t0, self.T, self.grid_size)
        self.solution = np.zeros([len(self.t), 13])
        #
        self.fitting_error_DF = 0.0
        self.fitting_error_DHF = 0.0
        self.r_zero_cond = False
        self.fitting_error_DF_cond = False
        self.fitting_error_DHF_cond = False
        self.stop_condition = False

    def update_conditions_search(self):
        r_zero_cond = (self.r_zero > 1)
        fitting_error_DF_cond = (self.fitting_error_DF < self.bound_error_FD)
        fitting_error_DHF_cond = (self.fitting_error_DHF
                                  < self.bound_error_FHD)
        stop_condition = (fitting_error_DF_cond and fitting_error_DHF_cond) \
                         and r_zero_cond

        self.stop_condition = stop_condition
        self.r_zero_cond = r_zero_cond
        self.fitting_error_DF_cond = fitting_error_DF_cond
        self.fitting_error_DHF_cond = fitting_error_DHF_cond
        return stop_condition

    def fitting_plot(self):
        t = self.t
        Y_m1_h = self.solution[:, 8]
        z = self.solution[:, 12]
        t_data_DF = self.frecuency_per_week_DF[3:, 0]
        t_data_DHF = self.frecuency_per_week_DHF[1:, 0]
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
        delte_index_t_Y = [0, 1, 3, 4, 6, 7, 8]
        t_Y_m1_h = np.delete(t_Y_m1_h, delte_index_t_Y)
        Y_m1_h_points = Y_m1_h[0: -1: offset]
        Y_m1_h_points = np.delete(Y_m1_h_points, delte_index_t_Y)

        #
        frecuency_per_week_DF = self.frecuency_per_week_DF[3:, 1]
        frecuency_per_week_DHF = self.frecuency_per_week_DHF[1:, 1]
        fitting_error_DF = \
            np.linalg.norm(frecuency_per_week_DF - z_points, ord=np.inf)
        self.fitting_error_DF = fitting_error_DF
        fitting_error_DHF = \
            np.linalg.norm(frecuency_per_week_DHF - Y_m1_h_points, ord=np.inf)
        self.fitting_error_DHF = fitting_error_DHF
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
        plt.savefig('./plots/fitting_DF_DHF.png')
