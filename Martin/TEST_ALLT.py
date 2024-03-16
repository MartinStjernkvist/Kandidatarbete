import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
from PLOT_Funktioner import Plot_functions
from FORM_Exempel import Ex
from FORM_8Performance import Perf
from FORM_General import Gen
from FORM_DesignTask2 import DT2
from FORM_Exempel import Ex
from KONSTANTER import *

T_11k = Gen().temperature(h_cruise)
a = Ex().a(gamma_air, R, T_11k)

P_1 = Gen().pressure(h_cruise)
T_1 = Gen().temperature(h_cruise)
a_1 = a
c_p_1 = C_p_air

"""---------------|Fan|---------------"""

P_2 = Ex().p_0(P_1, gamma_air, M_cruise)
T_2 = Ex().T_0(T_1, gamma_air, M_cruise)
gamma_air_2 = gamma_air
c_p_2 = c_p_1

r_hub_fan_nu = nu_fan * r_tip_fan[1]
ratio_hub_fan = (1 - (r_hub_fan_nu / r_hub_fan[1])) * 100  # procentuell skillnad

a_2 = Ex().a(gamma_air_2, R, T_2)  # speed of sound, station2
Ca_fan = M_ax_fan * a_2  # speed of air
u_fan = DT2().U(M_rel, a, Ca_fan)  # blade speed
omega_fan = u_fan / r_tip_fan[1]  # vinkelhastighet för fläkt
rpm_fan = DT2().RPM(omega_fan)

u_mid_fan = DT2().U_mid(r_tip_fan[1], r_hub_fan[1], omega_fan)  # medel av hastighet utöver fläkten
delta_h_fan = DT2().delta_H(psi_fan, u_mid_fan)
delta_T_fan = delta_h_fan / c_p_2

T_3 = T_2 + delta_T_fan  # temperaturen efter fläkten
FPR = Ex().stagnation_PR_fcn_TR(T_3, T_2, eta_fan, gamma_air_2)
P_3 = P_2 * FPR

MFP = Gen().MFP(M_2, gamma_air_2)
massflow = Gen().mdot(MFP, P_2, T_2, A_2)

"""---------------|LPC - 3 stage|---------------"""

gamma_air_3 = gamma_air
c_p_3 = C_p_air

omega_LPC = omega_fan  # pga samma axel, alltså samma vinkelhastighet

u_mid_LPC = []
for i in range(1, len(r_tip_LPC)):
    u_mid_LPC.append(DT2().U_mid(r_tip_LPC[i], r_hub_LPC[i], omega_LPC))
u_sum_LPC = sum(u ** 2 for u in u_mid_LPC)

delta_h_LPC = DT2().delta_H(psi_LPC, u_sum_LPC)  # designtask2, gånger 3 pga three stage LPC
delta_T_LPC = delta_h_LPC / c_p_3

T_4 = T_3 + delta_T_LPC
LPC_PR = Ex().stagnation_PR_fcn_TR(T_4, T_3, eta_LPC, gamma_air_3)
P_4 = P_3 * LPC_PR

"""---------------|HPC - 4 stage|---------------"""

gamma_air_4 = gamma_air
c_p_4 = C_p_air

a_4 = Ex().a(gamma_air_4, R, T_4)
Ca_HPC = M_ax_HPC * a_4

u_HPC = DT2().U(M_rel_HPC, a_4, Ca_HPC)
omega_HPC = u_HPC / r_tip_HPC
rmp_HPC = DT2().RPM(omega_HPC)

u_mid_HPC = []
for i in range(1, len(r_tip_HPC)):
    u_mid_HPC.append(DT2().U_mid(r_tip_HPC[i], r_hub_HPC[i], omega_HPC[i]))
u_sum_HPC = sum(u ** 2 for u in u_mid_HPC)

delta_h_HPC = DT2().delta_H(psi_HPC, u_sum_HPC)
delta_T_HPC = delta_h_HPC / c_p_4

T_5 = delta_T_HPC + T_4
HPC_PR = Ex().stagnation_PR_fcn_TR(T_5, T_4, eta_HPC, gamma_air_4)
P_5 = P_4 * HPC_PR

'''
lägg till eventuell kylning
'''

"""---------------|Combustor|---------------"""

gamma_air_5 = gamma_air
c_p_5 = C_p_air

c_p_mixed_5 = Perf().Cp_prod(T_5)
gamma_mixed_5 = Perf().gamma(c_p_mixed_5, R)  # Ändra senare efter rätt gamma mixed.

delta_h_combustion = c_p_mixed_5  # hur mkt bränsle

P_6 = P_5 * Combustor_ratio  # KOLLA SÅ RÄTT RATIO, tror rätt att trycket minskar?

"""---------------|HPT - 1 stage|---------------"""

omega_HPT = omega_HPC  # samma vinkelhastighet

gamma_mixed_6 = gamma_mixed_5
c_p_6 = C_p_air

u_HPT = r_tip_HPT * omega_LPC  # medel av hastighet utöver fläkten

u_sum_HPT = u_HPT ** 2

delta_h_HPT = DT2().delta_H(psi_HPT, u_sum_HPT)
delta_T_HPT = delta_h_HPT / c_p_6

T_7 = T_6 - delta_T_HPT

HPT_PR = Ex().stagnation_PR_fcn_TR(T_7, T_6, eta_HPT, gamma_mixed_6)
P_7 = P_6 * HPT_PR

"""---------------|LPT - 3 stage|---------------"""

c_p_7 = Perf().Cp_prod(T_7)
gamma_mixed_7 = Perf().gamma(c_p_7, R)
omega_LPT = omega_fan  # samma axel

u_mid_LPT = []
for i in range(1, len(r_tip_LPT)):
    u_mid_LPT.append(r_tip_LPT[i] * omega_LPT)

u_sum_LPT = sum(u ** 2 for u in u_mid_LPT)

delta_h_LPT = DT2().delta_H(psi_LPT, u_sum_LPT)
delta_T_LPT = delta_h_LPT / c_p_7

T_8 = T_7 - delta_T_LPT

LPT_PR = Ex().stagnation_PR_fcn_TR(T_8, T_7, eta_LPT, gamma_mixed_7)
P_8 = P_7 * LPT_PR

r_ratio_BPR = r_tip_fan[1] / r_tip_LPT[3]
BPR = ((np.log(1.9176 - (r_ratio_BPR * 1.25))) / (-0.2503)) - 0.6410

print(f'BPR = {BPR}\nMassflow = {massflow}')
