import numpy as np
import datetime
import yaml
from scipy import integrate
import matplotlib.pyplot as plt


class NumericOdeSolution:
    def __init__(self):
        self.Lambda_M = 5990.428571428572
        self.Lambda_S = 1.400000
        self.Lambda_m1 = 0.015000
        self.Lambda_m2 = 0.065000
        self.beta_M = 0.0305
        self.beta_H = 0.0329
        self.b = 2.095346
        self.mu_M = 0.09505714285714285
        self.mu_H = 0.000039
        self.alpha_c = 0.55340
        self.alpha_h = 0.55340  # 0.256219
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
        self.Y_m1_h0 = 400
        self.Y_m1_c0 = 400
        self.Y_m2_c0 = 0.0
        self.Rec_0 = 0.0
        self.z0 = 1.050000
        self.t0 = 0.0
        self.T = 250.000000
        self.grid_size = 2500000
        self.h = self.T / self.grid_size
        self.r_01 = 0.0
        self.r_02 = 0.0
        self.r_zero = 0
        #
        self.t = np.linspace(0, self.T, 100000)
        self.solution = np.zeros([len(self.t), 13])

    def f_rhs(self, x, t, Lambda_M, Lambda_S, Lambda_m1, Lambda_m2,
              beta_M, beta_H, b, mu_M, mu_H, alpha_c, alpha_h, sigma, p, 
              theta):
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
        S_m2 = x[7]
        Y_m1_h = x[8]
        Y_m1_c = x[9]
        Y_m2_c = x[10]
        R = x[11]
        z = x[12]
        # Load model parameters

        # Infection Forces
        N_H = I_s + I_1 + I_2 + S_m1 + S_m2 + Y_m1_h + Y_m1_c + Y_m2_c + R
        # N_M = M_s + M_I1 + M_I2
        c_M = (beta_M * b / N_H)
        c_H = (beta_H * b / N_H)
        A_I1 = c_M * I_1
        A_I2 = c_M * I_2
        A_Ym1_h = c_M * Y_m1_h
        A_Ym1_c = c_M * Y_m1_c
        A_Ym2_c = c_M * Y_m2_c
        B_M1 = c_H * M_I1
        B_M2 = c_H * M_I2

        # rhs of the ODE model
        dM_s = Lambda_M - (A_I1 + A_I2 + A_Ym1_c + A_Ym1_h + A_Ym2_c) * M_s \
               - mu_M * M_s
        dM_I1 = (A_I1 + A_Ym2_c) * M_s - mu_M * M_I1
        dM_I2 = (A_I2 + A_Ym1_h + A_Ym1_c) * M_s - mu_M * M_I2
        #
        dI_s = Lambda_S - (B_M1 + B_M2) * I_s - mu_H * I_s
        dI_1 = B_M1 * I_s - (alpha_c + mu_H) * I_1
        dI_2 = B_M2 * I_s - (alpha_c + mu_H) * I_2
        #
        dS_m1 = Lambda_m1 - sigma * B_M2 * S_m1 - mu_H * S_m1
        dS_m2 = Lambda_m2 - B_M1 * S_m2 - mu_H * S_m2
        #
        dY_m1_h = theta * sigma * B_M2 * S_m1 - (alpha_h + mu_H) * Y_m1_h
        dY_m1_c = (1.0 - theta) * sigma * B_M2 * S_m1 - (alpha_c + mu_H) * \
                                                       Y_m1_c
        dY_m2_c = B_M1 * S_m2 - (alpha_c + mu_H) * Y_m2_c
        #
        dR = alpha_c * (
            I_1 + I_2 + Y_m1_c + Y_m2_c) + alpha_h * Y_m1_h - mu_H * R
        dz = p * (dI_1 + dI_2 + dY_m1_c + dY_m2_c)
        dydt = np.array([dM_s, dM_I1, dM_I2,
                         dI_s, dI_1, dI_2,
                         dS_m1, dS_m2,
                         dY_m1_h, dY_m1_c, dY_m2_c,
                         dR, dz])
        dydt = dydt.astype('float64')
        return dydt

    def ode_int_solution(self):
        T = self.T
        t0 = self.t0
        t = np.linspace(t0, T, self.grid_size)
        y_0 = np.array(
            [self.M_s0, self.M_10, self.M_20,
             self.I_s0, self.I_10, self.I_20,
             self.S_m10, self.S_m20,
             self.Y_m1_h0, self.Y_m1_c0,
             self.Y_m2_c0, self.Rec_0,
             self.z0])
        Lambda_M = self.Lambda_M
        Lambda_S = self.Lambda_S
        Lambda_m1 = self.Lambda_m1
        Lambda_m2 = self.Lambda_m2
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
        #
        y = integrate.odeint(self.f_rhs, y_0, t, 
                             args=(Lambda_M, Lambda_S, Lambda_m1, Lambda_m2,
                                   beta_M, beta_H, b, mu_M, mu_H, alpha_c, 
                                   alpha_h, sigma, p, theta))
        self.solution = y
        self.t = t
        return y

    def parameters_sampling(self):
        # Mosquitoes constants
        Lambda_M = 20000
        Lambda_m1 = .015
        Lambda_m2 = 0.065
        beta_M = 0.05 * np.random.rand()
        b = abs(3 + np.random.randn())
        mu_M = np.random.randint(10, 22) ** (-1)

        # Human constants
        Lambda_S = 1.4
        beta_H = 0.05 * np.random.rand()
        mu_H = (365 * np.random.randint(70, 75)) ** (-1)

        # for recovering
        alpha_c = 0.215 + .05 * np.random.rand()
        alpha_h = .2 * alpha_c

        # Data fit parameters
        sigma = 5
        p = 0.05
        theta = .04 * np.random.rand()

        # Numerical parameters
        h = 0.1
        T = 200

        # Initial condition mosquitoes
        M_s0 = 120000
        M_10 = 20
        M_20 = 30
        #
        # Initial condition humans
        I_s0 = 35600
        I_10 = 1
        I_20 = 20
        S_m10 = 4400
        S_m20 = 0
        Y_m1_c0 = 0
        Y_m1_h0 = 0
        Y_m2_c0 = 0
        Rec_0 = 0
        z0 = p * (I_10 + I_20 + Y_m1_c0 + Y_m2_c0)

        # object parameters update
        self.Lambda_M = Lambda_M
        self.Lambda_S = Lambda_S
        self.Lambda_m1 = Lambda_m1
        self.Lambda_m2 = Lambda_m2
        self.beta_M = beta_M
        self.beta_H = beta_H
        self.b = b
        self.mu_M = mu_M
        self.mu_H = mu_H
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
        self.S_m10 = S_m10
        self.S_m20 = S_m20
        self.Y_m1_c0 = Y_m1_c0
        self.Y_m1_h0 = Y_m1_h0
        self.Y_m2_c0 = Y_m2_c0
        self.Rec_0 = Rec_0
        self.z0 = z0
        self.h = h
        self.T = T

        #
        new_parameters = [Lambda_M, Lambda_S, Lambda_m1, Lambda_m2, beta_M,
                          beta_H, b, mu_M, mu_H, alpha_c, alpha_h,
                          sigma, p, theta, M_s0, M_10, M_20, I_s0, I_10, I_20,
                          S_m10, S_m20, Y_m1_c0, Y_m1_h0, Y_m2_c0,
                          Rec_0, z0, h, T]
        new_parameters = np.array(new_parameters)
        return new_parameters

    def save_parameters(self,
                        file_name_prefix='./OutputParameters/parameters'):
        # load parameters
        Lambda_M = self.Lambda_M
        Lambda_S = self.Lambda_S
        Lambda_m1 = self.Lambda_m1
        Lambda_m2 = self.Lambda_m2
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
        S_m10 = self.S_m10
        S_m20 = self.S_m20
        Y_m1_c0 = self.Y_m1_c0
        Y_m1_h0 = self.Y_m1_h0
        Y_m2_c0 = self.Y_m2_c0
        Rec_0 = self.Rec_0
        z0 = self.z0
        h = self.h
        T = self.T
        parameters = {
            'Lambda_M': Lambda_M, 'Lambda_S': Lambda_S,
            'Lambda_m1': Lambda_m1, 'Lambda_m2': Lambda_m2,
            'beta_M': beta_M, 'beta_H': beta_H, 'b': b,
            'mu_M': mu_M, 'mu_H': mu_H, 'alpha_c': alpha_c,
            'alpha_h': alpha_h, 'sigma': sigma, 'p': p,
            'theta': theta, 'M_s0': M_s0, 'M_10': M_10,
            'M_20': M_20, 'I_s0': I_s0, 'I_10': I_10,
            'I_20': I_20, 'S_m10': S_m10, 'S_m20': S_m20,
            'Y_m1_c0': Y_m1_c0, 'Y_m1_h0': Y_m1_h0,
            'Y_m2_c0': Y_m2_c0, 'Rec_0': Rec_0, 'z0': z0, 'h': h,
            'T': T
            }

        str_time = str(datetime.datetime.now())
        file_name = file_name_prefix + str_time + '.yml'

        with open(file_name, 'w') as outfile:
            yaml.dump(parameters, outfile, default_flow_style=False)

    def compute_r_zero(self):
        # load parameters
        Lambda_M = self.Lambda_M
        Lambda_S = self.Lambda_S
        Lambda_m1 = self.Lambda_m1
        Lambda_m2 = self.Lambda_m2
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
        #
        N_M = Lambda_M / mu_M
        N_H = (Lambda_S + Lambda_m1 + Lambda_m2) / mu_H
        #
        pi_r = (beta_H * beta_M * b ** 2 * N_M) / (N_H ** 2 * mu_M * mu_H)
        r_01 = (pi_r / (alpha_c + mu_H)) * (Lambda_S + Lambda_m2)
        #
        r_02 = pi_r * (
            Lambda_S / (alpha_c + mu_H) + (theta * sigma * Lambda_m1) / (
                alpha_h + mu_H)
            + (1.0 - theta) * sigma * Lambda_m1 / (alpha_c + mu_H))

        self.r_01 = np.sqrt(r_01)
        self.r_02 = np.sqrt(r_02)
        r_zero = np.max(np.array([self.r_01, self.r_02]))
        self.r_zero = r_zero
        return np.sqrt(r_01), np.sqrt(r_02), r_zero

    def solution_plot(self):

        M_s = self.solution[:, 0]
        M_1 = self.solution[:, 1]
        M_2 = self.solution[:, 2]
        I_s = self.solution[:, 3]
        I_1 = self.solution[:, 4]
        I_2 = self.solution[:, 5]
        S_m1 = self.solution[:, 6]
        S_m2 = self.solution[:, 7]
        Y_m1_h = self.solution[:, 8]
        Y_m1_c = self.solution[:, 9]
        Y_m2_c = self.solution[:, 10]
        recovers = self.solution[:, 11]
        z = self.solution[:, 12]
        #
        t = self.t
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
        ax_array[2, 1].plot(t, S_m2)
        ax_array[2, 1].set_title(r'$S_{-2}$')
        ax_array[2, 2].plot(t, recovers)
        ax_array[2, 2].set_title(r'Recovered')

        #
        ax_array[3, 0].plot(t, Y_m1_h)
        ax_array[3, 0].set_title(r'$Y_{-1h}$')

        ax_array[3, 1].plot(t, Y_m1_c)
        ax_array[3, 1].set_title(r'$Y_{-1c}$')

        ax_array[3, 2].plot(t, Y_m2_c)
        ax_array[3, 2].set_title(r'$Y_{-2c}$')

        for i in np.arange(3):
            ax_array[3, i].set(xlabel='time')
        for j in np.arange(4):
            ax_array[j, 0].set(ylabel='Individuals')
        #
        plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
        plt.savefig('./plots/populations_grid.png')
        plt.close(f1)

        plt.figure(2)
        #
        plt.subplot(211)
        plt.plot(t, z,
                 ls='-',
                 color='blue',
                 alpha=0.4
                 )
        plt.xlabel(r'time(weeks)')
        plt.ylabel(r'$p * (I_1 + I_2 + Y_{-1c} + Y_{-2c})$')
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
        plt.savefig('./plots/DF_DHF.eps')
        plt.close(2)
        # plt.show()

    def load_parameters(self, file_name):
        with open(file_name, 'r') as f:
            parameter_data = yaml.load(f)
        # Set initial conditions

        self.Lambda_M = np.float64(parameter_data.get('Lambda_M'))
        self.Lambda_S = np.float64(parameter_data.get('Lambda_S'))
        self.Lambda_m1 = np.float64(parameter_data.get('Lambda_m1'))
        self.Lambda_m2 = np.float64(parameter_data.get('Lambda_m2'))
        self.beta_M = np.float64(parameter_data.get('beta_M'))
        self.beta_H = np.float64(parameter_data.get('beta_H'))
        self.b = np.float64(parameter_data.get('b'))
        self.mu_M = np.float64(parameter_data.get('mu_M'))
        self.mu_H = np.float64(parameter_data.get('mu_H'))
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
        self.S_m10 =np.float64(parameter_data.get('S_m10'))
        self.S_m20 = np.float64(parameter_data.get('S_m20'))
        self.Y_m1_h0 = np.float64(parameter_data.get('Y_m1_h0'))
        self.Y_m1_c0 = np.float64(parameter_data.get('Y_m1_c0'))
        self.Y_m2_c0 = np.float64(parameter_data.get('Y_m2_c0'))
        self.Rec_0 = np.float64(parameter_data.get('Rec_0'))
        self.z0 = np.float64(parameter_data.get('z0'))
        self.h = np.float64(parameter_data.get('h'))
        self.T = np.float64(parameter_data.get('T'))
