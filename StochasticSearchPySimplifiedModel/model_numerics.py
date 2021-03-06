import numpy as np
import datetime
import yaml
from scipy import integrate
import matplotlib.pyplot as plt


class NumericOdeSolution:
    def __init__(self):
        self.Lambda_M = 5990.428571428572
        self.beta_M = 0.0305
        self.beta_H = 0.0329
        self.b = 2.095346
        self.mu_M = 0.09505714285714285
        self.mu_H = 0.000039
        self.alpha_c = 0.55340
        self.alpha_h = 0.55340  # 0.256219
        self.kappa_1 = 0.7
        self.kappa_2 = 0.3
        self.sigma = 5.000000
        self.p = 0.050000
        self.theta = 0.450000
        self.M_s0 = 120000.000000
        self.M_10 = 200.000000
        self.M_20 = 300.000000
        self.I_s0 = 3600.000000
        self.I_e0 = 0.0
        self.I_10 = 100.000000
        self.I_20 = 200.000000
        self.S_m1_0 = 4400.000000
        self.Y_m1_0 = 400
        self.Rec_0 = 0.0
        self.z0 = .05 * (self.I_10 + self.I_20 + self.Y_m1_0)
        self.t0 = 0.0
        self.T = 250.000000
        self.grid_size = 2500000
        self.h = self.T / self.grid_size
        self.r_01 = 0.0
        self.r_02 = 0.0
        self.r_zero = 0
        self.N_H = 20000.0
        #
        self.t = np.linspace(0, self.T, 100000)
        self.solution = np.zeros([len(self.t), 13])

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

    def f_exposed_rhs(self, x, t,Lambda_M, beta_M, beta_H, b, mu_M, alpha_c,
                                   alpha_h, kappa_1, kappa_2, sigma, p, theta):
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

        # rhs of the ODE model
        dM_s = Lambda_M - (A_I1 + A_I2 + A_Ym1) * M_s \
               - mu_M * M_s
        dM_I1 = A_I1 * M_s - mu_M * M_I1
        dM_I2 = (A_I2 + A_Ym1) * M_s - mu_M * M_I2
        #
        dI_s = -(B_M1 + B_M2) * I_s
        dI_e = (B_M1 + B_M2) * I_s - (kappa_1 + kappa_2) * I_e
        dI_1 = kappa_1 * I_e - alpha_c * I_1
        dI_2 = kappa_2 * I_e - alpha_c * I_2
        #
        dS_m1 = - sigma * B_M2 * S_m1
        #
        dY_m1 = sigma * B_M2 * S_m1 - alpha_c * Y_m1
        #
        dR = alpha_c * (I_1 + I_2 + Y_m1)
        dz = 0.05 * (dI_1 + dI_2 + (1.0 - theta) * dY_m1)
        # dz = p * (kappa_1 + kappa_2) * I_e
        dY_m1_h = theta * sigma * B_M2 * S_m1 - alpha_h * Y_m1_h
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
        kappa_1 = self.kappa_1
        kappa_2 = self.kappa_2
        sigma = self.sigma
        p = self.p
        theta = self.theta
        #
        y = integrate.odeint(self.f_exposed_rhs, y_0, t,
                             args=(Lambda_M, beta_M, beta_H, b, mu_M, alpha_c,
                                   alpha_h, kappa_1, kappa_2, sigma, p, theta))
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
        theta = .04 + 0.5 * np.random.rand()

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

