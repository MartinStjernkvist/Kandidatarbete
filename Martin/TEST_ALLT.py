import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
from FORM_Exempelupg import Ex
from FORM_8Performance import Perf
from FORM_General import Gen
from FORM_DesignTask2 import DT2
from FORM_Exempelupg import Ex
from KONSTANTER import *

"""---------------|Innan motor|---------------"""

T_cruise = Gen().temperature(h_cruise)
a = Ex().a(gamma_air, R, T_cruise)

# P_1 = Gen().pressure(h_cruise)
P_1 = 7171
T_1 = Gen().temperature(h_cruise)
a_1 = a
C_p_1 = C_p_air

"""---------------|Fan|---------------"""

P_2 = Ex().p_0(P_1, gamma_air, M_cruise)
T_2 = Ex().T_0(T_1, gamma_air, M_cruise)

gamma_air_2 = gamma_air
C_p_2 = C_p_1

# Jämförelse av hubben genom hub tip ratio och mätta värden
r_hub_fan_nu = nu_fan * r_tip_fan[1]
ratio_hub_fan = (1 - (r_hub_fan_nu / r_hub_fan[1])) * 100  # procentuell skillnad

a_2 = Ex().a(gamma_air_2, R, T_2)
c_air_fan = M_ax_fan * a_2  # speed of air

U_fan = DT2().U(M_rel, a, c_air_fan)  # blade speed
omega_fan = U_fan / r_tip_fan[1]  # vinkelhastighet för fläkt
RPM_fan = Gen().RPM(omega_fan)

# Tar inte hänsyn till konisk form på hubben
U_mid_fan = DT2().U_mid(r_tip_fan[1], r_hub_fan[1], omega_fan)  # medelhastighet över fläkt
sum_U_squared_fan = U_mid_fan ** 2

# Steglast och flödesfaktor
Delta_H_fan = DT2().delta_H(psi_fan, sum_U_squared_fan)
Delta_T_fan = Delta_H_fan / C_p_2

T_3 = T_2 + Delta_T_fan  # temperaturen efter fläkten
FPR = Ex().stagnation_PR_fcn_TR(T_3, T_2, eta_fan, gamma_air_2)
P_3 = P_2 * FPR

MFP = Gen().MFP(M_2, gamma_air_2)
mass_flow = Gen().mass_flow(MFP, P_2, T_2, A_2)

"""---------------|LPC - 3 stage|---------------"""

gamma_air_3 = gamma_air
C_p_3 = C_p_air

omega_LPC = omega_fan  # pga samma axel, alltså samma vinkelhastighet

U_mid_LPC = []
for i in range(1, len(r_tip_LPC)):
    U_mid_LPC.append(DT2().U_mid(r_tip_LPC[i], r_hub_LPC[i], omega_LPC))

sum_U_squared_LPC = sum(u ** 2 for u in U_mid_LPC)

Delta_H_LPC = DT2().delta_H(psi_LPC, sum_U_squared_LPC)  # designtask2, gånger 3 pga three stage LPC
Delta_T_LPC = Delta_H_LPC / C_p_3

T_4 = T_3 + Delta_T_LPC
LPC_PR = Ex().stagnation_PR_fcn_TR(T_4, T_3, eta_LPC, gamma_air_3)
P_4 = P_3 * LPC_PR

"""---------------|HPC - 4 stage|---------------"""

gamma_air_4 = gamma_air
C_p_4 = C_p_air

a_4 = Ex().a(gamma_air_4, R, T_4)
c_air_HPC = M_ax_HPC * a_4

U_HPC = DT2().U(M_rel_HPC, a_4, c_air_HPC)
omega_HPC = U_HPC / r_tip_HPC
RPM_HPC = Gen().RPM(omega_HPC)

U_mid_HPC = []
for i in range(1, len(r_tip_HPC)):
    U_mid_HPC.append(DT2().U_mid(r_tip_HPC[i], r_hub_HPC[i], omega_HPC[i]))

sum_U_squared_HPC = sum(u ** 2 for u in U_mid_HPC)

Delta_H_HPC = DT2().delta_H(psi_HPC, sum_U_squared_HPC)
Delta_T_HPC = Delta_H_HPC / C_p_4

T_5 = T_4 + Delta_T_HPC
HPC_PR = Ex().stagnation_PR_fcn_TR(T_5, T_4, eta_HPC, gamma_air_4)

cooling = 0.8
P_5 = P_4 * HPC_PR * cooling

# Tillägg kylning, ca 20% Ta enbart ut kylning efter HPC, sida 104-105 Mattingly.

'''
lägg till eventuell kylning
'''

"""---------------|Combustor|---------------"""

# gamma_air_5 = gamma_air
# C_p_5 = C_p_air
T_6 = T_turbine

m_c = 0.05 + 0.00015 * (T_6 - 1200) # eqn (19) 8Performance

# NYTT GAMMA OCH CP FÖR GAS OCH AIR
C_p_mix_5 = Perf().C_p_prod(T_5)
gamma_mix_5 = Perf().gamma(C_p_mix_5, R)

Delta_H_combustion = C_p_mix_5  # hur mkt bränsle

P_6 = P_5 * combustor_PR  # KOLLA SÅ RÄTT RATIO, tror rätt att trycket minskar?

# massflow_mixed = massflow_ny + massflow_fuel

"""---------------|HPT - 1 stage|---------------"""

gamma_mixed_6 = gamma_mix_5
C_p_mix_6 = Perf().C_p_prod(T_6)

omega_HPT = omega_HPC

U_HPT = r_tip_HPT * omega_LPC

sum_U_squared_HPT = U_HPT ** 2

Delta_H_HPT = DT2().delta_H(psi_HPT, sum_U_squared_HPT)
Delta_T_HPT = Delta_H_HPT / C_p_mix_6

T_7 = T_6 + Delta_T_HPT

HPT_PR = Ex().stagnation_PR_fcn_TR(T_7, T_6, eta_HPT, gamma_mixed_6)
P_7 = P_6 * HPT_PR

"""---------------|LPT - 3 stage|---------------"""

C_p_mix_7 = Perf().C_p_prod(T_7)
gamma_mix_7 = Perf().gamma(C_p_mix_7, R)

omega_LPT = omega_fan  # samma axel

U_mid_LPT = []
for i in range(1, len(r_tip_LPT)):
    U_mid_LPT.append(r_tip_LPT[i] * omega_LPT)

sum_U_squared_LPT = sum(u ** 2 for u in U_mid_LPT)

Delta_H_LPT = DT2().delta_H(psi_LPT, sum_U_squared_LPT)
Delta_T_LPT = Delta_H_LPT / C_p_mix_7

T_8 = T_7 + Delta_T_LPT

LPT_PR = Ex().stagnation_PR_fcn_TR(T_8, T_7, eta_LPT, gamma_mix_7)
P_8 = P_7 * LPT_PR

r_ratio_BPR = r_tip_fan[1] / r_tip_LPT[3]
BPR = ((np.log(1.9176 - (r_ratio_BPR * 1.25))) / (-0.2503)) - 0.6410

# print(f'BPR = {BPR}\nMassflow = {mass_flow}')

parameter_list = [
    P_2,
    T_2,
    r_hub_fan_nu,
    ratio_hub_fan,
    a_2,
    c_air_fan,
    U_fan,
    omega_fan,
    RPM_fan,
    U_mid_fan,
    sum_U_squared_fan,
    Delta_H_fan,
    Delta_T_fan,
    T_3,
    FPR,
    P_3,
    MFP,
    mass_flow,
    C_p_3,
    omega_LPC,
    U_mid_LPC,
    sum_U_squared_LPC,
    Delta_H_LPC,
    Delta_T_LPC,
    T_4,
    LPC_PR,
    P_4,
    gamma_air_4,
    C_p_4,
    a_4,
    c_air_HPC,
    U_HPC,
    omega_HPC,
    RPM_HPC,
    U_mid_HPC,
    sum_U_squared_HPC,
    Delta_H_HPC,
    Delta_T_HPC,
    T_5,
    HPC_PR,
    P_5,
    gamma_mix_5,
    C_p_mix_5,
    C_p_mix_5,
    gamma_mix_5,
    Delta_H_combustion,
    P_6,
    omega_HPT,
    gamma_mixed_6,
    C_p_mix_6,
    U_HPT,
    sum_U_squared_HPT,
    Delta_H_HPT,
    Delta_T_HPT,
    T_6,
    T_7,
    HPT_PR,
    P_7,
    C_p_mix_7,
    gamma_mix_7,
    omega_LPT,
    U_mid_LPT,
    sum_U_squared_LPT,
    Delta_H_LPT,
    Delta_T_LPT,
    T_8,
    LPT_PR,
    P_8,
    r_ratio_BPR,
    BPR
]

name_list = [
    'OK: P_2',
    'OK: T_2',
    'r_hub_fan_nu',
    '*** NOT OK ***: ratio_hub_fan',
    'a_2',
    'c_air_fan',
    'U_fan',
    'omega_fan',
    'RPM_fan',
    'U_mid_fan',
    'sum_U_squared_fan',
    'Delta_H_fan',
    'Delta_T_fan',
    'T_3',
    'FPR',
    'P_3',
    'MFP',
    'mass_flow',
    'C_p_3',
    'omega_LPC',
    'U_mid_LPC',
    'sum_U_squared_LPC',
    'Delta_H_LPC',
    'Delta_T_LPC',
    'T_4',
    'LPC_PR',
    'P_4',
    'gamma_air_4',
    'C_p_4',
    'a_4',
    'c_air_HPC',
    'U_HPC',
    'omega_HPC',
    'RPM_HPC',
    'U_mid_HPC',
    'sum_U_squared_HPC',
    'Delta_h_HPC',
    'Delta_T_HPC',
    'OK: T_5',
    'HPC_PR',
    'P_5',
    'gamma_air_5',
    'C_p_5',
    'C_p_mixed_5',
    'gamma_mixed_5',
    'Delta_H_combustion',
    'P_6',
    'omega_HPT',
    'gamma_mixed_6',
    'C_p_6',
    'U_HPT',
    'sum_U_squared_HPT',
    'Delta_H_HPT',
    'Delta_T_HPT',
    'T_6',
    'T_7',
    'HPT_PR',
    'P_7',
    'C_p_7',
    'gamma_mixed_7',
    'omega_LPT',
    'U_mid_LPT',
    'sum_U_squared_LPT',
    'Delta_H_LPT',
    'Delta_T_LPT',
    'T_8',
    'LPT_PR',
    'P_8',
    'r_ratio_BPR',
    'BPR'
]
values = [(name_list[item], parameter_list[item]) for item in range(len(parameter_list))]
for i in range(len(values)):
    print(values[i])
