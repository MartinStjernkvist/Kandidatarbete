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

    def calc_theta_0break(self, T_t4max, T_tSLS):
        return T_t4max / T_tSLS


class AppendixD(CommonFunctionalityAppD):
    def __init__(self, T_t4max, T_tSLS, gamma_c, eta_c, eta_m, c_pt, c_pc, T_0, M_0):
        self.T_t4max = T_t4max
        self.T_tSLS = T_tSLS
        self.gamma_c = gamma_c
        self.eta_c = eta_c
        self.eta_m = eta_m
        self.c_pt = c_pt
        self.c_pc = c_pc
        self.M_0 = M_0
        self.T_0 = T_0
        self.theta_0_break = self.calc_theta_0break(T_t4max, T_tSLS)

    def theta_0(self):
        """
        Equation (D.1)
        Same as (2.52a), but with gamma_c
        """
        theta_tau_r = (self.T_0 / T_std) * (1 + ((self.gamma_c -1) / 2) * self.M_0 ** 2)
        return theta_tau_r


    def CompressorRatio(self):
        return (1 + )

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
