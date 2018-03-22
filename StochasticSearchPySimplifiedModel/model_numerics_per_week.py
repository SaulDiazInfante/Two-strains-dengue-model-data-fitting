import numpy as np
from scipy import integrate
from model_numerics import NumericOdeSolution
import matplotlib.pyplot as plt


class NumericOdeSolutionPerWeek(NumericOdeSolution):

    def __init__(self):
        self.Lambda_M = 7 * 20000
        self.beta_M = 0.029149
        self.beta_H = 0.030855
        self.b = 7 * 2.095346
        self.mu_M = 7 * 0.076923
        self.mu_H = 7 * 0.000039
        self.alpha_c = 7 * 0.256219
        self.alpha_h = 7 * 0.256219
        self.kappa_1 = .7
        self.kappa_2 = 0.3
        self.sigma = 5
        self.p = 0.05
        self.theta = 0.019
        self.M_s0 = 120000.000000
        self.M_10 = 20.000000
        self.M_20 = 30.000000
        self.I_s0 = 35600.000000
        self.I_e0 = 30.0
        self.I_10 = 10.000000
        self.I_20 = 20.000000
        self.S_m1_0 = 4400.000000
        self.Y_m1_0 = 25.0
        self.Rec_0 = 0.0
        self.z0 = self.p * (self.I_10 + self.I_20 + self.Y_m1_0)
        self.T = 52.000000
        self.t0 = 25.0
        self.grid_size = np.int(self.T - self.t0) * 10000
        self.h = self.T / self.grid_size
        self.r_01 = 0.0
        self.r_02 = 0.0
        self.r_zero = 0
        self.N_H = 20000.0
        #
        self.t = np.linspace(0, self.T, 100000)
        self.solution = np.zeros([len(self.t), 13])

    def parameters_sampling_per_week(self, flag_deterministic=False):

        if flag_deterministic:
            # Mosquitoes constants
            Lambda_M = 7 * 20000
            beta_M = 0.05  # + 0.45 * np.random.rand()
            b = 7 * .5
            mu_M = 7 * 0.02

            # Human constants
            beta_H = 0.05
            # for recovering
            alpha_c = 0.083  # + 0.17 * np.random.rand()) * 7
            alpha_h = alpha_c  # 0.125 * 7
            sigma = 2.5  # + 4.0 * np.random.rand()
            kappa_1 = 0.83
            kappa_2 = 0.83
            p = 0.05
            theta = 0.01  # + 0.4 * np.random.rand()

            # Initial condition mosquitoes
            M_s0 = 100
            M_10 = 0
            M_20 = 0

            # Initial condition humans

            # self.N_H = 20000 + 10000 * np.random.rand()
            I_s0 = self.N_H
            I_10 = 0
            I_20 = 0
            S_m1_0 = 0
            #
            Y_m1_0 = 0
            Rec_0 = 0.0
            z0 = p * (I_10 + I_20 + Y_m1_0)
        else:
            Lambda_M = 7 * np.abs(6000 + 2000 * np.random.randn())
            beta_M = 0.05 * np.random.rand()
            b = np.abs(2 + np.random.randn()) * 7
            mu_M = 7 * (0.033 + 0.067 * np.random.rand())

            # Human constants
            beta_H = 0.05 * np.random.rand()
            # for recovering
            alpha_c = (0.0556 + .00444 * np.random.rand()) * 7
            alpha_h = (0.125 + 0.125 * np.random.rand()) * 7
            kappa_1 = 7 * (0.05 * (0.17 - 0.083) + np.abs(np.random.randn()))
            kappa_2 = 7 * (0.05 * (0.17 - 0.083) + np.abs(np.random.randn()))
            sigma = 0.5 + 4.5 * np.random.rand()
            p = 0.05
            theta = 0.01 + 0.025 * np.random.rand()

            # Initial condition mosquitoes
            M_s0 = 0.70 * (Lambda_M / mu_M)
            M_10 = 0.0 # 0.05 * (Lambda_M / mu_M)
            M_20 = 0.0 # 0.05 * (Lambda_M / mu_M)

            # Initial condition humans
            self.N_H = 30000 + 10000 * np.random.rand()
            I_s0 = 0.9 * self.N_H
            I_10 = 0.025 * self.N_H
            I_20 = 0.025 * self.N_H
            S_m1_0 = 0.045 * self.N_H
            #
            Y_m1_0 = 0.005 * self.N_H
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
        self.kappa_1 = kappa_1
        self.kappa_2 = kappa_2
        self.sigma = sigma
        self.p = p
        self.theta = theta
        #
        self.M_s0 = M_s0
        self.M_10 = M_10
        self.M_20 = M_20
        #
        self.I_s0 = .8 * self.N_H
        self.I_10 = .15 * self.N_H
        self.I_20 = .05 * self.N_H
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
                          beta_H, b, mu_M, alpha_c, alpha_h, kappa_1, kappa_2,
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
        kappa_1 = self.kappa_1
        kappa_2 = self.kappa_2
        #
        y = integrate.odeint(self.f_exposed_rhs, y_0, t,
                             args=(Lambda_M, beta_M, beta_H, b, mu_M, alpha_c,
                                   alpha_h, sigma, p, theta, kappa_1, kappa_2))
        self.solution = y
        self.t = t
        return y
