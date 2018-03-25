import numpy as np
import datetime
import yaml
from scipy import integrate
import matplotlib.pyplot as plt


class StochasticSearch():
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
        self.kappa = 0.17
        self.rho_1 = .8
        self.rho_2 = .1
        self.M_s0 = 120000.000000
        self.M_10 = 20.000000
        self.M_20 = 30.000000
        self.I_s0 = 35600.000000
        self.I_e0 = 30.0
        self.I_10 = 1.000000
        self.I_20 = 20.000000
        self.S_m1_0 = 4400.000000
        self.Y_m1_0 = 100.0
        self.Rec_0 = 0.0
        self.z0 = 1.050000
        self.t0 = 25  # 25.0
        self.T = 53
        self.grid_size = np.int(self.T - self.t0) * 10000
        self.h = np.float64(self.T) / np.float64(self.grid_size)
        self.r_01 = 0.0
        self.r_02 = 0.0
        self.r_zero = 0
        self.N_H = 20000.0
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
        stop_condition = (fitting_error_DF_cond or fitting_error_DHF_cond) \
                        and r_zero_cond

        self.stop_condition = stop_condition
        self.r_zero_cond = r_zero_cond
        self.fitting_error_DF_cond = fitting_error_DF_cond
        self.fitting_error_DHF_cond = fitting_error_DHF_cond
        return stop_condition

    def fitting_plot(self):
        t = self.t
        # z = self.solution[:, 10]
        # Y_m1_h = self.solution[:, 11]
        Y_m1_h = self.solution[:, 10]
        z = self.solution[:, 9]
        #
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

    def fitting_error(self):
        t = self.t
        # z = self.solution[:, 10]
        # Y_m1_h = self.solution[:, 11]
        Y_m1_h = self.solution[:, 10]
        z = self.solution[:, 9]
        #
        t_data_DF = self.frecuency_per_week_DF[3:, 0]
        t_data_DHF = self.frecuency_per_week_DHF[1:, 0]
        offset = 10000
        #
        t_z = t[0: -1: offset]
        t_z = np.round(t_z)
        t_z = t_z.astype(int)
        z_points = z[0:-1: offset]
        #
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
        #
        frecuency_per_week_DF = self.frecuency_per_week_DF[3:, 1]
        frecuency_per_week_DHF = self.frecuency_per_week_DHF[1:, 1]
        fitting_error_DF = \
            np.linalg.norm(frecuency_per_week_DF - z_points, ord=np.inf) \
            / np.linalg.norm(frecuency_per_week_DF, ord=np.inf)
        self.fitting_error_DF = fitting_error_DF
        fitting_error_DHF = \
            np.linalg.norm(frecuency_per_week_DHF - Y_m1_h_points,
                           ord=np.inf) / \
            np.linalg.norm(frecuency_per_week_DHF, ord=np.inf)
        self.fitting_error_DHF = fitting_error_DHF

    def fitting_error_exposed(self):

        t = self.t
        z = self.solution[:, 10]
        Y_m1_h = self.solution[:, 11]
        #
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
            np.linalg.norm(frecuency_per_week_DF - z_points, ord=np.inf) / \
            np.linalg.norm(frecuency_per_week_DF, ord=np.inf)
        self.fitting_error_DF = fitting_error_DF
        fitting_error_DHF = \
            np.linalg.norm(frecuency_per_week_DHF \
                           - Y_m1_h_points, ord=np.inf) / \
            np.linalg.norm(frecuency_per_week_DHF, ord=np.inf)
        self.fitting_error_DHF = fitting_error_DHF
        #
        #

    def fitting_plot_exposed(self):

        t = self.t
        z = self.solution[:, 10]
        Y_m1_h = self.solution[:, 11]
        #
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

    def f_rhs(self, x, t, Lambda_M, beta_M, beta_H, b, mu_M,
              alpha_c, alpha_h, sigma, p, theta):
        """

        :type alpha_h: float64
        """
        M_s = x[0]
        M_I1 = x[1]
        M_I2 = x[2]
        I_s = x[3]
        I_1 = x[4]
        I_2 = x[5]
        S_m1 = x[6]
        Y_m1 = x[7]
        R = x[8]
        z = x[9]
        Y_m1_h = x[10]
        # Load model parameters

        # Infection Forces
        N_H = self.N_H
        # N_M = M_s + M_I1 + M_I2
        c_M = (beta_M * b / N_H)
        c_H = (beta_H * b / N_H)
        A_I1 = c_M * I_1
        A_I2 = c_M * I_2
        A_Ym1 = c_M * Y_m1
        B_M1 = c_H * M_I1
        B_M2 = c_H * M_I2

        # rhs of the ODE model
        dM_s = Lambda_M - (A_I1 + A_I2 + A_Ym1) * M_s \
               - mu_M * M_s
        dM_I1 = A_I1 * M_s - mu_M * M_I1
        dM_I2 = (A_I2 + A_Ym1) * M_s - mu_M * M_I2
        #
        dI_s = -(B_M1 + B_M2) * I_s
        dI_1 = B_M1 * I_s - alpha_c * I_1
        dI_2 = B_M2 * I_s - alpha_c * I_2
        #
        dS_m1 = - sigma * B_M2 * S_m1
        #
        dY_m1 = sigma * B_M2 * S_m1 - (alpha_c + alpha_h) * Y_m1
        #
        dR = alpha_c * (I_1 + I_2 + Y_m1) + alpha_h * Y_m1
        dz = p * (dI_1 + dI_2 + (1.0 - theta) * dY_m1)
        dY_m1_h = theta * dY_m1
        dydt = np.array([dM_s, dM_I1, dM_I2,
                         dI_s, dI_1, dI_2,
                         dS_m1, dY_m1, dR, dz, dY_m1_h])
        dydt = dydt.astype('float64')
        return dydt

    def f_exposed_rhs(self, x, t, Lambda_M, beta_M, beta_H, b, mu_M, alpha_c,
                      alpha_h, kappa, sigma, p, theta):
        """

        :param x:
        :param t:
        :param Lambda_M:
        :param beta_M:
        :param beta_H:
        :param b:
        :param mu_M:
        :param alpha_c:
        :param alpha_h:
        :param sigma:
        :param p:
        :param theta:
        :param kappa_1:
        :param kappa_2:
        :type alpha_h: float64
    """
        M_s = x[0]
        M_I1 = x[1]
        M_I2 = x[2]
        I_s = x[3]
        I_e = x[4]
        I_1 = x[5]
        I_2 = x[6]
        S_m1 = x[7]
        Y_m1 = x[8]
        R = x[9]
        z = x[10]
        Y_m1_h = x[11]

        # Load model parameters

        # Infection Forces
        N_H = self.N_H
        # N_M = M_s + M_I1 + M_I2
        c_M = (beta_M * b / N_H)
        c_H = (beta_H * b / N_H)
        A_I1 = c_M * I_1
        A_I2 = c_M * I_2
        A_Ym1 = c_M * Y_m1
        B_M1 = c_H * M_I1
        B_M2 = c_H * M_I2
        rho_1 = self.rho_1
        rho_2 = self.rho_2
        # rhs of the ODE model
        dM_s = Lambda_M - (A_I1 + A_I2 + A_Ym1) * M_s \
               - mu_M * M_s
        dM_I1 = A_I1 * M_s - mu_M * M_I1
        dM_I2 = (A_I2 + A_Ym1) * M_s - mu_M * M_I2
        #
        dI_s = -(B_M1 + B_M2) * I_s
        dI_e = (B_M1 + B_M2) * I_s + sigma * B_M2 * S_m1 - kappa * I_e
        dI_1 = rho_1 * kappa * I_e - alpha_c * I_1
        dI_2 = rho_2 * kappa * I_e - alpha_c * I_2
        #
        dS_m1 = - sigma * B_M2 * S_m1
        #
        dY_m1 = (1.0 - (rho_1 + rho_2)) * kappa * I_e - alpha_c * Y_m1
        #
        dR = alpha_c * (I_1 + I_2 + Y_m1)
        dz = 0.05 * (dI_1 + dI_2 + (1.0 - theta) * dY_m1)
        dY_m1_h = theta * dY_m1
        dydt = np.array([dM_s, dM_I1, dM_I2,
                         dI_s, dI_e, dI_1, dI_2,
                         dS_m1, dY_m1, dR, dz, dY_m1_h])
        dydt = dydt.astype('float64')
        return dydt

    def ode_int_solution(self):
        T = self.T
        t0 = self.t0
        t = np.linspace(t0, T, self.grid_size)
        y_0 = np.array(
            [self.M_s0, self.M_10, self.M_20,
             self.I_s0, self.I_10, self.I_20,
             self.S_m1_0, self.Y_m1_0,
             self.Rec_0, self.z0, self.theta * self.Y_m1_0])
        Lambda_M = self.Lambda_M
        beta_M = self.beta_M
        beta_H = self.beta_H
        b = self.b
        mu_M = self.mu_M
        alpha_c = self.alpha_c
        alpha_h = self.alpha_h
        sigma = self.sigma
        p = self.p
        theta = self.theta
        #
        y = integrate.odeint(self.f_rhs, y_0, t,
                             args=(Lambda_M, beta_M, beta_H, b, mu_M, alpha_c,
                                   alpha_h, sigma, p, theta))
        self.solution = y
        self.t = t
        return y

    def ode_int_solution_exposed(self):
        T = self.T
        t0 = self.t0
        t = np.linspace(t0, T, self.grid_size)
        y_0 = np.array(
            [self.M_s0, self.M_10, self.M_20,
             self.I_s0, self.I_e0, self.I_10, self.I_20,
             self.S_m1_0, self.Y_m1_0,
             self.Rec_0, self.z0, self.theta * self.Y_m1_0])
        Lambda_M = self.Lambda_M
        beta_M = self.beta_M
        beta_H = self.beta_H
        b = self.b
        mu_M = self.mu_M
        alpha_c = self.alpha_c
        alpha_h = self.alpha_h
        kappa = self.kappa
        sigma = self.sigma
        p = self.p
        theta = self.theta
        #
        y = integrate.odeint(self.f_exposed_rhs, y_0, t,
                             args=(Lambda_M, beta_M, beta_H, b, mu_M, alpha_c,
                                   alpha_h, kappa, sigma, p, theta))
        self.solution = y
        self.t = t
        return y

    def parameters_sampling(self):
        # Mosquitoes constants
        Lambda_M = 40000 + 30000 * np.random.rand()
        beta_M = 0.1 + np.random.rand()
        b = abs(4 + np.random.randn())
        mu_M = np.random.randint(10, 22) ** (-1)

        beta_H = 0.25 + np.random.rand()
        # for recovering
        alpha_c = 0.215 + .05 * np.random.rand()
        alpha_h = .4 * alpha_c

        # Data fit parameters
        sigma = 1.5 + 3.5 * np.random.rand()
        p = 0.05
        theta = .04 + 0.05 * np.random.rand()

        # Numerical parameters
        h = 0.1
        T = 200

        # Initial condition mosquitoes
        M_s0 = 120000
        M_10 = 1200
        M_20 = 600
        #
        # Initial condition humans
        I_s0 = 35600
        I_10 = 10
        I_20 = 20
        S_m1_0 = 4400
        Y_m1_0 = 500 + 1000 * np.random.rand()
        Rec_0 = 0
        z0 = p * (I_10 + I_20 + Y_m1_0)

        # object parameters update
        self.Lambda_M = Lambda_M
        self.beta_M = beta_M
        self.beta_H = beta_H
        self.b = b
        self.mu_M = mu_M
        self.alpha_c = alpha_c
        self.alpha_h = alpha_h
        self.sigma = sigma
        self.p = p
        self.theta = theta
        #
        self.M_s0 = M_s0
        self.M_10 = M_10
        self.M_20 = M_20
        self.I_s0 = I_s0
        self.I_10 = I_10
        self.I_20 = I_20
        self.S_m1_0 = S_m1_0
        self.Y_m1_0 = Y_m1_0
        self.Rec_0 = Rec_0
        self.z0 = z0
        self.h = h
        self.T = T
        #
        #
        new_parameters = [Lambda_M, beta_M, beta_H, b, mu_M, alpha_c,
                          alpha_h, sigma, p, theta, M_s0, M_10, M_20, I_s0,
                          I_10, I_20, S_m1_0, Y_m1_0,
                          Rec_0, z0, h, T]
        new_parameters = np.array(new_parameters)
        return new_parameters

    def save_parameters(self,
                        file_name_prefix='./OutputParameters/parameters'):
        # load parameters
        Lambda_M = self.Lambda_M
        beta_M = self.beta_M
        beta_H = self.beta_H
        b = self.b
        mu_M = self.mu_M
        mu_H = self.mu_H
        alpha_c = self.alpha_c
        alpha_h = self.alpha_h
        sigma = self.sigma
        p = self.p
        theta = self.theta
        M_s0 = self.M_s0
        M_10 = self.M_10
        M_20 = self.M_20
        I_s0 = self.I_s0
        I_10 = self.I_10
        I_20 = self.I_20
        S_m1_0 = self.S_m1_0
        Y_m1_0 = self.Y_m1_0
        Rec_0 = self.Rec_0
        z0 = self.z0
        h = self.h
        T = self.T
        parameters = {
            'Lambda_M': Lambda_M,
            'beta_M': beta_M, 'beta_H': beta_H, 'b': b,
            'mu_M': mu_M, 'alpha_c': alpha_c,
            'alpha_h': alpha_h, 'sigma': sigma, 'p': p,
            'theta': theta, 'M_s0': M_s0, 'M_10': M_10,
            'M_20': M_20, 'I_s0': I_s0, 'I_10': I_10,
            'I_20': I_20, 'S_m1_0': S_m1_0, 'Y_m1_0': Y_m1_0,
            'Rec_0': Rec_0, 'z0': z0, 'h': h,
            'T': T
            }

        str_time = str(datetime.datetime.now())
        file_name = file_name_prefix + str_time + '.yml'
        with open(file_name, 'w') as outfile:
            yaml.dump(parameters, outfile, default_flow_style=False)

    def compute_r_zero(self):
        # load parameters
        Lambda_M = self.Lambda_M
        beta_M = self.beta_M
        beta_H = self.beta_H
        b = self.b
        mu_M = self.mu_M
        alpha_c = self.alpha_c
        alpha_h = self.alpha_h
        sigma = self.sigma
        I_s0 = self.I_s0
        S_m1_0 = self.S_m1_0
        # p = self.p
        # theta = self.theta
        #
        N_M = Lambda_M / mu_M
        N_H = self.N_H
        #
        pi_r = (beta_H * beta_M * b ** 2 * Lambda_M) / (N_H ** 2 * mu_M ** 2)
        r_01 = pi_r * I_s0 * alpha_c ** (-1)
        #
        r_02 = pi_r * sigma * S_m1_0 * (alpha_c + alpha_h) ** (-1)

        self.r_01 = r_01
        self.r_02 = r_02
        r_zero = np.sqrt(self.r_01 + self.r_02)
        self.r_zero = r_zero
        return np.sqrt(r_01), np.sqrt(r_02), r_zero

    def solution_plot(self):
        """
        M_s = self.solution[:, 0]
        M_1 = self.solution[:, 1]
        M_2 = self.solution[:, 2]
        I_s = self.solution[:, 3]
        I_e = self.solution[:, 4]
        I_1 = self.solution[:, 5]
        I_2 = self.solution[:, 6]
        S_m1 = self.solution[:, 7]
        Y_m1 = self.solution[:, 8]
        recovers = self.solution[:, 9]
        z = self.solution[:, 10]
        Y_m1_h = self.solution[:, 11]
        """
        M_s = self.solution[:, 0]
        M_1 = self.solution[:, 1]
        M_2 = self.solution[:, 2]
        I_s = self.solution[:, 3]
        I_1 = self.solution[:, 4]
        I_2 = self.solution[:, 5]
        S_m1 = self.solution[:, 6]
        Y_m1 = self.solution[:, 7]
        recovers = self.solution[:, 8]
        z = self.solution[:, 9]
        Y_m1_h = self.solution[:, 10]
        #
        t = self.t
        N_H = I_s + I_1 + I_2 + S_m1 + Y_m1 + recovers
        #
        f1, ax_array = plt.subplots(4, 3, sharex=True)
        #
        ax_array[0, 0].plot(t, M_s)
        ax_array[0, 0].set_title(r'$M_s$')

        ax_array[0, 1].plot(t, M_1)
        ax_array[0, 1].set_title(r'$M_1$')

        ax_array[0, 2].plot(t, M_2)
        ax_array[0, 2].set_title(r'$M_2$')
        #
        #   Infected humans first time
        ax_array[1, 0].plot(t, I_s)
        ax_array[1, 0].set_title(r'$I_s$')

        ax_array[1, 1].plot(t, I_1)
        ax_array[1, 1].set_title(r'$I_1$')

        ax_array[1, 2].plot(t, I_2)
        ax_array[1, 2].set_title(r'$I_2$')
        ##
        ##
        ax_array[2, 0].plot(t, S_m1)
        ax_array[2, 0].set_title(r'$S_{-1}$')
        ax_array[2, 1].plot(t, Y_m1)
        ax_array[2, 1].set_title(r'$Y_{-1}$')
        ax_array[2, 2].plot(t, recovers)
        ax_array[2, 2].set_title(r'Recovered')
        #
        #
        ax_array[3, 0].plot(t, N_H)
        ax_array[3, 0].set_title(r'$N_H$')

        ax_array[3, 1].plot(t, Y_m1_h)
        ax_array[3, 1].set_title(r'$Y_{-1_h}$')
        #
        ax_array[3, 2].plot(t, z)
        ax_array[3, 2].set_title(r'$z$')
        #
        for i in np.arange(3):
            ax_array[3, i].set(xlabel='time')
        for j in np.arange(4):
            ax_array[j, 0].set(ylabel='Individuals')
        #
        plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
        plt.savefig('./plots/populations_grid.png')
        plt.close(f1)
        #
        #
        plt.figure(2)
        #
        plt.subplot(211)
        plt.plot(t, z,
                 ls='-',
                 color='blue',
                 alpha=0.4
                 )
        plt.xlabel(r'time(weeks)')
        plt.ylabel(r'$p * (I_1 + I_2 + Y_{-1})$')
        #
        plt.subplot(212)
        plt.plot(t, Y_m1_h,
                 ls='-',
                 color='red',
                 alpha=0.4
                 )
        plt.xlabel(r'time(weeks)')
        plt.ylabel(r'$Y_{-1h}$')
        plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
        # plt.show()
        plt.savefig('./plots/DF_DHF.png')
        plt.close(2)

    def solution_plot_exposed(self):
        M_s = self.solution[:, 0]
        M_1 = self.solution[:, 1]
        M_2 = self.solution[:, 2]
        I_s = self.solution[:, 3]
        I_e = self.solution[:, 4]
        I_1 = self.solution[:, 5]
        I_2 = self.solution[:, 6]
        S_m1 = self.solution[:, 7]
        Y_m1 = self.solution[:, 8]
        recovers = self.solution[:, 9]
        z = self.solution[:, 10]
        Y_m1_h = self.solution[:, 11]
        #
        t = self.t
        N_H = I_s + +I_e + I_1 + I_2 + S_m1 + Y_m1 + recovers
        #
        f1, ax_array = plt.subplots(4, 3, sharex=True)
        #
        ax_array[0, 0].plot(t, M_s)
        ax_array[0, 0].set_title(r'$M_s$')

        ax_array[0, 1].plot(t, M_1)
        ax_array[0, 1].set_title(r'$M_1$')

        ax_array[0, 2].plot(t, M_2)
        ax_array[0, 2].set_title(r'$M_2$')
        #
        #   Infected humans first time
        ax_array[1, 0].plot(t, I_s)
        ax_array[1, 0].set_title(r'$I_s$')

        ax_array[1, 1].plot(t, I_1)
        ax_array[1, 1].set_title(r'$I_1$')

        ax_array[1, 2].plot(t, I_2)
        ax_array[1, 2].set_title(r'$I_2$')
        ##
        ##
        ax_array[2, 0].plot(t, S_m1)
        ax_array[2, 0].set_title(r'$S_{-1}$')
        ax_array[2, 1].plot(t, Y_m1)
        ax_array[2, 1].set_title(r'$Y_{-1}$')
        ax_array[2, 2].plot(t, recovers)
        ax_array[2, 2].set_title(r'Recovered')
        #
        #
        ax_array[3, 0].plot(t, N_H)
        ax_array[3, 0].set_title(r'$N_H$')
        #
        ax_array[3, 1].plot(t, Y_m1_h)
        ax_array[3, 1].set_title(r'$Y_{-1_h}$')
        #
        ax_array[3, 2].plot(t, z)
        ax_array[3, 2].set_title(r'$z$')
        #
        for i in np.arange(3):
            ax_array[3, i].set(xlabel='time')
        for j in np.arange(4):
            ax_array[j, 0].set(ylabel='Individuals')
        #
        plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
        plt.savefig('./plots/populations_grid.png')
        plt.close(f1)
        #
        #
        plt.figure(2)
        #
        plt.subplot(211)
        plt.plot(t, z,
                 ls='-',
                 color='blue',
                 alpha=0.4
                 )
        plt.xlabel(r'time(weeks)')
        plt.ylabel(r'$p * (I_1 + I_2 + Y_{-1})$')
        #
        plt.subplot(212)
        plt.plot(t, Y_m1_h,
                 ls='-',
                 color='red',
                 alpha=0.4
                 )
        plt.xlabel(r'time(weeks)')
        plt.ylabel(r'$Y_{-1h}$')
        plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
        # plt.show()
        plt.savefig('./plots/DF_DHF.png')
        plt.close(2)

    def load_parameters(self, file_name):
        with open(file_name, 'r') as f:
            parameter_data = yaml.load(f)
        # Set initial conditions
        #
        self.Lambda_M = np.float64(parameter_data.get('Lambda_M'))
        self.beta_M = np.float64(parameter_data.get('beta_M'))
        self.beta_H = np.float64(parameter_data.get('beta_H'))
        self.b = np.float64(parameter_data.get('b'))
        self.mu_M = np.float64(parameter_data.get('mu_M'))
        self.alpha_c = np.float64(parameter_data.get('alpha_c'))
        self.alpha_h = np.float64(parameter_data.get('alpha_h'))
        self.sigma = np.float64(parameter_data.get('sigma'))
        self.p = np.float64(parameter_data.get('p'))
        self.theta = np.float64(parameter_data.get('theta'))
        self.M_s0 = np.float64(parameter_data.get('M_s0'))
        self.M_10 = np.float64(parameter_data.get('M_10'))
        self.M_20 = np.float64(parameter_data.get('M_20'))
        self.I_s0 = np.float64(parameter_data.get('I_s0'))
        self.I_10 = np.float64(parameter_data.get('I_10'))
        self.I_20 = np.float64(parameter_data.get('I_20'))
        self.S_m1_0 = np.float64(parameter_data.get('S_m1_0'))
        self.Y_m1_0 = np.float64(parameter_data.get('Y_m1_0'))
        self.Rec_0 = np.float64(parameter_data.get('Rec_0'))
        self.z0 = np.float64(parameter_data.get('z0'))
        self.h = np.float64(parameter_data.get('h'))
        self.T = np.float64(parameter_data.get('T'))

    def parameters_sampling_per_week(self, flag_deterministic=False):

        if flag_deterministic:
            # Mosquitoes constants
            Lambda_M = 7 * 20000
            beta_M = 0.05
            b = 7 * .5
            mu_M = 7 * 0.02

            # Human constants
            beta_H = 0.05
            # for recovering
            alpha_c = 0.083
            alpha_h = alpha_c
            sigma = 2.5

            p = 0.05
            theta = 0.01
            #
            # Initial condition mosquitoes
            M_s0 = 100
            M_10 = 0
            M_20 = 0
            #
            # Initial condition humans
            # self.N_H = 20000 + 10000 * np.random.rand()
            I_s0 = self.N_H
            I_10 = 0
            I_20 = 0
            S_m1_0 = 0
            #
            #
            Y_m1_0 = 0
            Rec_0 = 0.0
            z0 = p * (I_10 + I_20 + Y_m1_0)
        else:
            Lambda_M = 7 * np.abs(6000 + 2000 * np.random.randn())
            beta_M = 0.05 * np.random.rand()
            b = np.abs(2 + np.random.randn()) * 7
            mu_M = 7 * (0.033 + 0.067 * np.random.rand())
            #
            # Human constants
            beta_H = 0.05 * np.random.rand()
            # for recovering
            alpha_c = (0.0556 + .0444 * np.random.rand()) * 7
            alpha_h = (0.125 + 0.125 * np.random.rand()) * 7
            sigma = 0.5 + 4.5 * np.random.rand()
            p = 0.05
            theta = 0.1 + 0.75 * np.random.rand()
            rho_1 = .6 + .2 * np.random.rand()
            rho_2 = .8 - rho_1
            # Initial condition mosquitoes
            p1 = 0.9 + .1 * np.random.rand()
            pj = np.random.rand(2)
            pj_hat = .9 * (1.0 - p1) / pj.sum() * pj

            M_s0 = p1 * (Lambda_M / mu_M)
            M_10 = pj_hat[0] * (Lambda_M / mu_M)
            M_20 = pj_hat[1] * (Lambda_M / mu_M)

            # Initial condition humans
            self.N_H = 37000 + 5000 * np.random.rand()
            # partition
            p1 = 0.8 + .15 * np.random.rand()
            pj = 1 - p1

            I_10 = 1.0
            I_20 = 1.0
            Y_m1_0 = 1.0
            #
            #
            I_s0 = p1 * (self.N_H - (I_10 + I_20 + Y_m1_0))
            S_m1_0 = pj * self.N_H
            #
            #
            Rec_0 = 0.0
            z0 = p * (I_10 + I_20 + Y_m1_0 * (1 - theta))

        # Numerical parameters
        T = self.T
        h = np.float64(self.T) / np.float64(self.grid_size)
        #
        # object parameters update
        self.Lambda_M = Lambda_M
        self.beta_M = beta_M
        self.beta_H = beta_H
        self.b = b
        self.mu_M = mu_M
        self.alpha_c = alpha_c
        self.alpha_h = alpha_h
        self.sigma = sigma
        self.p = p
        self.theta = theta
        self.rho_1 = rho_1
        self.rho_2 = rho_2
        #
        self.M_s0 = M_s0
        self.M_10 = M_10
        self.M_20 = M_20
        #
        self.I_s0 = I_s0
        self.I_10 = I_10
        self.I_20 = I_20
        #
        self.S_m1_0 = S_m1_0
        #
        self.Y_m1_0 = Y_m1_0
        self.Rec_0 = Rec_0
        self.z0 = z0
        #
        self.h = h
        self.T = T
        #
        new_parameters = [Lambda_M, beta_M,
                          beta_H, b, mu_M, alpha_c, alpha_h,
                          sigma, p, theta, M_s0, M_10, M_20, I_s0, I_10, I_20,
                          S_m1_0, Y_m1_0, Rec_0, z0, h, T]
        new_parameters = np.array(new_parameters)
        return new_parameters

    def ode_int_solution_per_week(self):
        T = self.T
        t0 = self.t0
        t = np.linspace(t0, T, self.grid_size)
        y_0 = np.array(
            [self.M_s0, self.M_10, self.M_20,
             self.I_s0, self.I_10, self.I_20,
             self.S_m1_0, self.Y_m1_0,
             self.Rec_0, self.z0, self.theta * self.Y_m1_0])
        Lambda_M = self.Lambda_M
        beta_M = self.beta_M
        beta_H = self.beta_H
        b = self.b
        mu_M = self.mu_M
        alpha_c = self.alpha_c
        alpha_h = self.alpha_h
        sigma = self.sigma
        p = self.p
        theta = self.theta
        #
        y = integrate.odeint(self.f_rhs, y_0, t,
                             args=(Lambda_M, beta_M, beta_H, b, mu_M,
                                   alpha_c, alpha_h, sigma, p, theta))
        self.solution = y
        self.t = t
        return y

    def ode_int_solution_exposed_peer_week(self):
        T = self.T
        t0 = self.t0
        t = np.linspace(t0, T, self.grid_size)
        y_0 = np.array(
            [self.M_s0, self.M_10, self.M_20,
             self.I_s0, self.I_e0, self.I_10, self.I_20,
             self.S_m1_0, self.Y_m1_0,
             self.Rec_0, self.z0, self.theta * self.Y_m1_0])
        Lambda_M = self.Lambda_M
        beta_M = self.beta_M
        beta_H = self.beta_H
        b = self.b
        mu_M = self.mu_M
        alpha_c = self.alpha_c
        alpha_h = self.alpha_h
        sigma = self.sigma
        p = self.p
        theta = self.theta
        kappa = self.kappa
        #
        y = integrate.odeint(self.f_exposed_rhs, y_0, t,
                             args=(Lambda_M, beta_M, beta_H, b, mu_M, alpha_c,
                                   alpha_h, sigma, p, theta, kappa))
        self.solution = y
        self.t = t
        return y
