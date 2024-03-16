import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
from PLOT_Funktioner import Plot_functions
from FORM_Exempel import Ex
from FORM_8Performance import Perf
from FORM_Mattingly import Gen
from FORM_Exempel import Ex
from KONSTANTER import *

P_1 = Ex().p_0(Gen().pressure(h_cruise), gamma_air, M_cruise)
T_1 = Ex().T_0(Gen().temperature(h_cruise), gamma_air, M_cruise)
C_P1 = Perf().Cp_air(T_1)

pi_LPC = np.sqrt(pi_C)
pi_HPC = pi_LPC
tau_LPC = Perf().tau(pi_LPC, eta_C)
T_2 = tau_LPC * P_1
P_2 = pi_LPC * P_1
C_P2 = Perf().Cp_air(T_2)
W_LPC = Perf().W(C_P1, C_P2, T_1, T_2)

T_3 = Perf().T_3(T_1, T_2, epsilon_IC)
P_3 = pi_IC * P_2
C_P3 = Perf().Cp_air(T_3)
tau_HPC = Perf().tau(pi_HPC, eta_C)

T_4 = tau_HPC * T_3
P_4 = pi_HPC * P_3
C_P4 = Perf().Cp_air(T_4)
W_HPC = Perf().W(C_P3, C_P4, T_3, T_4)

P_10 = P_1
P_9 = P_10 / pi_reg
P_5 = pi_burner * P_4
P_6 = pi_burner * P_5

pi_HPT = Perf().pi_HPT(P_9, P_6, pi_burner)
pi_LPT = pi_HPT

tau_HPT = Perf().tau_HPT(pi_LPT, eta_turbine)
tau_LPT = tau_HPT

T_6 = T_turbine
C_P6 = Perf().Cp_prod(T_6)
T_8 = T_turbine
C_P8 = Perf().Cp_prod(T_8)
T_7 = tau_HPT * T_turbine
C_P7 = Perf().Cp_prod(T_7)
T_9 = T_7
C_P9 = Perf().Cp_prod(T_9)
T_5 = Perf().T_5(T_4, T_9, epsilon_reg)
C_P5 = Perf().Cp_air(T_5)

m_C = 0.05 + 0.00015 * (T_turbine - 1200)

f_1 = Perf().f_1(m_C, T_5, T_6, C_P5, C_P6, Q_r, eta_burner)
f_2 = Perf().f_2(m_C, T_7, T_8, C_P7, C_P8, Q_r, eta_burner, f_1)

T_7C = Perf().T_7C(m_C, C_P7, T_7, C_P4, T_4, eta_burner, Q_r, C_P6, T_6)
T_9C = Perf().T_9C(m_C, f_1, f_2, C_P9, T_9, C_P4, T_4)
T_5C = Perf().T_5(T_4, T_9C, epsilon_reg)

f_1C = Perf().f_1(m_C, T_5C, T_6, C_P5, C_P6, Q_r, eta_burner)
f_2C = Perf().f_2(m_C, T_7C, T_8, C_P7, C_P8, Q_r, eta_burner, f_1C)

W_LPT = Perf().W_LPT(eta_mech, m_C, f_1C, f_2, C_P8, T_8, C_P9, T_9C)
W_HPT = Perf().W_HPT(eta_mech, m_C, f_1C, C_P6, T_6, C_P7, T_7C)
W_t = W_HPT + W_LPT

W_net = Perf().W_net(W_HPT, W_LPT, W_HPC, W_LPC)
SFC = Perf().SFC(W_net, f_1C, f_2C)
eta = Perf().eta(W_net, f_1C, f_2C)

print(eta)