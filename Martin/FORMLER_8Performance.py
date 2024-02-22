import numpy as np
from KONSTANTER import *


class Perf:
    def tau(self, pi, eta_C, gamma=gamma):
        return 1 + (pi ** ((gamma - 1) / gamma) - 1) / eta_C

    def W(self, C_P_start, C_P_end, T_start, T_end):
        return C_P_end * T_end - C_P_start * T_start

    def T_3(self, T_1, T_2, epsilon_IC):
        return T_2 - epsilon_IC * (T_2 - T_1)

    def tau_HPC(self, pi_HPC, eta_C, gamma=gamma):
        return self.tau(pi_HPC, eta_C, gamma=gamma)

    def tau_HPT(self, pi_LPT, eta_t, gamma_t=gamma_t):
        return 1 - eta_t * (1 - pi_LPT ** ((gamma_t - 1) / gamma_t))

    def T_5(self, T_4, T_9, epsilon_reg):
        return T_4 + epsilon_reg * (T_9 - T_4)

    def f_1(self, m_C, T_5, T_6, C_P_5, C_P_6, Q_r, eta_b):
        return (1 - 2 * m_C) * (C_P_6 * T_6 - C_P_5 * T_5) / (eta_b * Q_r - C_P_6 * T_6)

    def f_2(self, m_C, T_7, T_8, C_P_7, C_P_8, Q_r, eta_b, f_1):
        return (1 - 2 * m_C) * (C_P_8 * T_8 - C_P_7 * T_7) / (eta_b * Q_r - C_P_8 * T_8)
