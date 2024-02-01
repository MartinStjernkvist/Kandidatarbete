import numpy as np


class CommonFunctionalityAppD:

    def calc_theta_0break(self, T_t4max, T_tSLS):
        return T_t4max / T_tSLS


class AppendixD(CommonFunctionalityAppD):
    def __init__(self, T_t4max, T_tSLS, gamma_c, eta_c, eta_m, c_pt, c_pc):
        self.T_t4max = T_t4max
        self.T_tSLS = T_tSLS
        self.gamma_c = gamma_c
        self.eta_c = eta_c
        self.eta_m = eta_m
        self.c_pt = c_pt
        self.c_pc = c_pc
        self.theta_0break = self.calc_theta_0break(T_t4max, T_tSLS)

    def theta_0(self):

    def CompressorRatio(self):
        return (1 + )

    def ThrottleRatio(self):
        """
        Equation (D.6)
        """
        theta_0break = self.T_t4max / self.T_tSLS
        return theta_0break

    def MachBreak(self):
        """
        Equation (D.7)
        """
        M_0break = np.sqrt((2 / (self.gamma_c - 1)) * (self.theta_0break - 1))
        return M_0break
