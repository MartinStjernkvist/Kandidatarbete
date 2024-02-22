import numpy as np
from KONSTANTER import *


class Perf:

    def tau(self, pi, eta_C, gamma=gamma_a):
        return 1 + (pi ** ((gamma - 1) / gamma) - 1) / eta_C

    def W(self, C_P_start, C_P_end, T_start, T_end):
        return C_P_end * T_end - C_P_start * T_start

    def T_3(self, T_1, T_2, epsilon_IC):
        return T_2 - epsilon_IC * (T_2 - T_1)

    def pi_HPT(self, P_9, P_6, pi_b):
        return np.sqrt(P_9 / (P_6 * pi_b))

    def tau_HPC(self, pi_HPC, eta_C, gamma=gamma_a):
        return self.tau(pi_HPC, eta_C, gamma=gamma)

    def tau_HPT(self, pi_LPT, eta_t, gamma_t=gamma_t):
        return 1 - eta_t * (1 - pi_LPT ** ((gamma_t - 1) / gamma_t))

    def T_5(self, T_4, T_9, epsilon_reg):
        return T_4 + epsilon_reg * (T_9 - T_4)

    def f_1(self, m_C, T_5, T_6, C_P5, C_P6, Q_r, eta_b):
        return (1 - 2 * m_C) * (C_P6 * T_6 - C_P5 * T_5) / (eta_b * Q_r - C_P6 * T_6)

    def f_2(self, m_C, T_7, T_8, C_P7, C_P8, Q_r, eta_b, f_1):
        return (1 - m_C + f_1) * (C_P8 * T_8 - C_P7 * T_7) / (eta_b * Q_r - C_P8 * T_8)

    def T_7C(self, m_C, C_P7, T_7, C_P4, T_4, eta_b, Q_r, C_P6, T_6):
        return ((1 - 2 * m_C) * C_P7 * T_7 + m_C * C_P4 * T_4) / (eta_b * Q_r - C_P6 * T_6)

    def T_9C(self, m_C, f_1, f_2, C_P9, T_9, C_P4, T_4):
        return ((1 - m_C + f_1 + f_2) * C_P9 * T_9 + m_C * C_P4 * T_4) / ((1 + f_1 + f_2) * C_P9)

    def W_LPT(self, eta_mech, m_C, f_1C, f_2, C_P8, T_8, C_P9, T_9C):
        return eta_mech * (1 - m_C + f_1C + f_2) * (C_P8 * T_8 - C_P9 * T_9C)

    def W_HPT(self, eta_mech, m_C, f_1C, C_P6, T_6, C_P7, T_7C):
        return eta_mech * (1 - 2 * m_C + f_1C) * (C_P6 * T_6 - C_P7 * T_7C)

    def W_net(self, W_HPT, W_LPT, W_HPC, W_LPC):
        return W_HPT + W_LPT - W_HPC - W_LPC

    def SFC(self, W_net, f_1C, f_2C):
        return (f_1C + f_2C) / W_net

    def eta(self, W_net, f_1C, f_2C):
        return W_net / (f_1C + f_2C)

    def Cp_air(self, T):
        return A_0_air + A_1_air * T + A_2_air * T ** 2 + A_3_air * T ** 3 + A_4_air * T ** 4 + A_5_air * T ** 5 + A_6_air * T ** 6 + A_7_air * T ** 7

    def Cp_prod(self, T):
        return A_0_prod + A_1_prod * T + A_2_prod * T ** 2 + A_3_prod * T ** 3 + A_4_prod * T ** 4 + A_5_prod * T ** 5 + A_6_prod * T ** 6 + A_7_prod * T ** 7

    def Cp_mix(self, T, f):
        return (self.Cp_air(T) + f * self.Cp_prod(T)) / (1 + f)