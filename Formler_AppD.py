import numpy as np

"""
Constants (keep track)
"""
g_0 = 9.81          #
rho_SL = 1.225      #
t_R = 3             # runtime on takeoff
s_TO = 500          #
T_std = 288.15      #
P_std = 101325      #
gamma = 1.4         #
R = 287             # J/kg*K
placeholder = 10 ** 20  # use as placeholder when creating instances of classes,
                        # not used in calculations of the instance

class CommonFunctionalityAppD:

    def calc_theta_0_break(self, T_t_4_max, T_t_SLS):
        return T_t_4_max / T_t_SLS


class AppendixD(CommonFunctionalityAppD):
    def __init__(self, T_t_4_max, T_t_SLS, gamma_c, eta_c, eta_m, beta, c_pt, c_pc, T_0, M_0):
        self.T_t_4_max = T_t_4_max
        self.T_t_SLS = T_t_SLS
        self.gamma_c = gamma_c
        self.eta_c = eta_c
        self.eta_m = eta_m
        self.c_pt = c_pt
        self.c_pc = c_pc
        self.M_0 = M_0
        self.T_0 = T_0
        self.beta = beta
        self.theta_0_break = self.calc_theta_0break(T_t_4_max, T_t_SLS)

    def theta_0(self):
        """
        Equation (D.1)
        Same as (2.52a), but with gamma_c
        """
        theta_tau_r = (self.T_0 / T_std) * (1 + ((self.gamma_c -1) / 2) * self.M_0 ** 2)
        return theta_tau_r


    def CompressorTotalTemperatureRatioToYield(self):
        """
        Equation (D.2)

        Depends only on: turbine total temperature ratio tau_t, throttle setting: T_t_4
        """
        tau_c = 1 + self.eta_m * (1 - self.beta) * (1 + f)
        return tau_c

    def ThrottleRatio(self):
        """
        Equation (D.6)
        """
        theta_0_break = self.T_t4max / self.T_tSLS
        return theta_0_break

    def MachBreak(self):
        """
        Equation (D.7)
        """
        M_0break = np.sqrt((2 / (self.gamma_c - 1)) * (self.theta_0_break - 1))
        return M_0break
