import numpy as np
from FORMLER_8Performance import Perf
from KONSTANTER import *









f_1 = Perf().f_1(m_C, T_5, T_6, C_P_5, C_P_6, Q_r, eta_b)
f_2 = Perf().f_2(m_C, T_7, T_8, C_P_7, C_P_8, Q_r, eta_b, f_1)
W_LPT = eta_mech * (1 - m_C + f_1 + f_2) * (C_P_8 * T_8 - C_P_9 * T_9_C)
W_HPT = eta_mech * (1 - 2 * m_C + f_1) * (C_P_6 * T_6 - C_P_7 * T_7_C)

W_t = W_HPT + W_LPT
W_net = W_HPT -
SFC = (f_1 + f_2) / W_net
eta = W_net /(f_1 + f_2)