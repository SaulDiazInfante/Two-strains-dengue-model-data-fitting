import numpy as np
from scipy import integrate
from model_numerics import NumericOdeSolution
import matplotlib.pyplot as plt

class NumericOdeSolutionPerWeek(NumericOdeSolution):

    def __init__(self):
        self.Lambda_M = 7 * 20000
        self.Lambda_S = 7 * 1.4
        self.Lambda_m1 = 7 * 0.015000
        self.Lambda_m2 = 7 * 0.065000
        self.beta_M = 0.029149
        self.beta_H = 0.030855
        self.b = 7 * 2.095346
        self.mu_M = 7 * 0.076923
        self.mu_H = 7 * 0.000039
        self.alpha_c = 7 * 0.256219
        self.alpha_h = 7 * 0.256219
        self.sigma = 5
        self.p = 0.05
        self.theta = 0.019
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
        self.z0 = self.p * (self.I_10 + self.I_20 + self.Y_m1_c0 +
                            self.Y_m2_c0)
        self.T = 52.000000
        self.t0 = 25.0
        self.grid_size = np.int(self.T - self.t0) * 10000
        self.h = self.T / self.grid_size
        self.r_01 = 0.0
        self.r_02 = 0.0
        self.r_zero = 0
        #
        self.t = np.linspace(0, self.T, 100000)
        self.solution = np.zeros([len(self.t), 13])

    def parameters_sampling_per_week(self):
        # Mosquitoes constants
        Lambda_M = 7 * np.abs(41933 + 0*2000 * np.random.randn())
        Lambda_m1 = 0.015 * 7
        Lambda_m2 = 0.065 * 7
        beta_M = 0.05 * np.random.rand()
        b = np.abs(2 + np.random.randn()) * 7
        mu_M = 7 * (0.033 + 0.067 * np.random.rand())
        # Human constants
        Lambda_S = 1.4 * 7
        beta_H = 0.05 * np.random.rand()

        mu_H = (365 * np.random.randint(70, 75)) ** (-1) * 7

        # for recovering
        alpha_c = (0.0556 + .00444 * np.random.rand()) * 7
        alpha_h = (0.125 + 0.125 * np.random.rand()) * 7

        # Data fit parameters
        sigma = 2.0 + 2.8 * np.random.rand()
        p = 0.05
        theta = 0.01 + 0.09 * np.random.rand()

        # Numerical parameters
        T = self.T
        h = np.float64(self.T) / np.float64(self.grid_size)

        # Initial condition mosquitoes
        M_s0 = 120000
        M_10 = 20
        M_20 = 30
        #
        # Initial condition humans
        I_s0 = 35600.0
        I_10 = 1.0
        I_20 = 20.0
        S_m10 = 4400.0
        S_m20 = 0.0
        #
        Y_m1_h0 = 0.0  # 10 + np.random.randint(1, 10)
        Y_m1_c0 = 0.0  # 5 + np.random.randint(1, 10)
        Y_m2_c0 = 0.0  # 10 + np.random.randint(1, 10)
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
        #
        self.I_s0 = I_s0
        self.I_10 = I_10
        self.I_20 = I_20
        #
        self.S_m10 = S_m10
        self.S_m20 = S_m20
        #
        self.Y_m1_h0 = Y_m1_h0
        self.Y_m1_c0 = Y_m1_c0
        self.Y_m2_c0 = Y_m2_c0
        self.Rec_0 = Rec_0
        self.z0 = z0
        #
        self.h = h
        self.T = T
        #
        new_parameters = [Lambda_M, Lambda_S, Lambda_m1, Lambda_m2, beta_M,
                          beta_H, b, mu_M, mu_H, alpha_c, alpha_h,
                          sigma, p, theta, M_s0, M_10, M_20, I_s0, I_10, I_20,
                          S_m10, S_m20, Y_m1_h0, Y_m1_c0, Y_m2_c0,
                          Rec_0, z0, h, T]
        new_parameters = np.array(new_parameters)
        return new_parameters

    def ode_int_solution_per_week(self):
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
