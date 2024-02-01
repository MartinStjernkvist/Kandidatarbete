import numpy as np


class ThrottleRatio():
    def __init__(self, T_t4max, T_tSLS, gamma_c):
        self.T_t4max = T_t4max
        self.T_tSLS = T_tSLS
        self.gamma_c = gamma_c
        self.theta_0break = self.T_t4max / self.T_tSLS

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
