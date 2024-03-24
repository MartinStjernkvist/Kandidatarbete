import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
from FORM_Exempelupg import Ex
from FORM_8Performance import Perf
from FORM_General import Gen
from FORM_DesignTask2 import DT2
from FORM_Exempelupg import Ex
from KONSTANTER import *
from scipy.optimize import fsolve
"""
Lite syntax

P = tryck
T = temperatur
a = ljudhastighet
C_p = värmekapacitet (konstant tryck)
C_v = värmekapacitet (konstant temperatur)
gamma = specifik värme ratio
omega = vinkelhastighet
U = bladhastighet
c = lufthastighet
MFP = mass flow parameter
psi = steglast
Delta_H = entalpiändring
Delta_T = temperaturändring
"""


"""---------------|Innan motor|---------------"""

T_cruise = Gen().temperature(h_cruise)
a = Ex().a(gamma_air, R, T_cruise)

# P_1 = Gen().pressure(h_cruise)
P_0 = 7171
T_0 = Gen().temperature(h_cruise)
a_1 = a
C_p_1 = C_p_air

"""---------------|Fan|---------------"""

P_t_2 = Ex().p_0(P_0, gamma_air, M_cruise)

T_t_2 = Ex().T_0(T_0, gamma_air, M_cruise)
T_2 = T_t_2 / Ex().stagnation_TR(gamma_air, M_2)

C_p_2 = Perf().C_p_air(T_t_2)  # Stagnations temperatur används för att man kan tänka sig att gasen bromsas till c = 0, värms upp och sedan höjs hastigheten igen
gamma_air_2 = Perf().gamma(C_p_2, R)

MFP = Gen().MFP(M_2, gamma_air_2)
mass_flow = Gen().mass_flow(MFP, P_t_2, T_2, A_2)

# Jämförelse av hubben genom hub tip ratio och mätta värden
r_hub_fan_nu = nu_fan * r_tip_fan[1]
ratio_hub_fan = (1 - (r_hub_fan_nu / r_hub_fan[1])) * 100  # procentuell skillnad

a_2 = Ex().a(gamma_air_2, R, T_2)
c_air_fan = M_ax_fan * a_2

U_fan = DT2().U(M_rel_fan, a, c_air_fan)
omega_fan = U_fan / r_tip_fan[1]
RPM_fan = Gen().RPM(omega_fan)

# Tar inte hänsyn till konisk form på hubben
U_mid_fan = DT2().U_mid(r_tip_fan[1], r_hub_fan[1], omega_fan)
sum_U_squared_fan = U_mid_fan ** 2

Delta_H_fan = DT2().delta_H(psi_fan, sum_U_squared_fan)
Delta_T_fan = Delta_H_fan / C_p_2

T_t_21 = T_t_2 + Delta_T_fan  # temperaturen efter fläkten
FPR = Ex().stagnation_PR_fcn_TR(T_t_21, T_t_2, eta_fan, gamma_air_2)
P_t_21 = P_t_2 * FPR

"""---------------|LPC - 3 stage|---------------"""

C_p_21 = Perf().C_p_air(T_t_21)
gamma_air_21 = Perf().gamma(C_p_21, R)

A_21 = Gen().A(r_tip_LPC[1], r_hub_LPC[1])

omega_LPC = omega_fan  # pga samma axel

U_mid_LPC = []
for i in range(1, len(r_tip_LPC)):
    U_mid_LPC.append(DT2().U_mid(r_tip_LPC[i], r_hub_LPC[i], omega_LPC))

sum_U_squared_LPC = sum(u ** 2 for u in U_mid_LPC)

Delta_H_LPC = DT2().delta_H(psi_LPC, sum_U_squared_LPC)  # gånger 3 pga three stage LPC
Delta_T_LPC = Delta_H_LPC / C_p_21

T_t_25 = T_t_21 + Delta_T_LPC
LPC_PR = Ex().stagnation_PR_fcn_TR(T_t_25, T_t_21, eta_LPC, gamma_air_21)
P_t_25 = P_t_21 * LPC_PR



"""---------------|HPC - 4 stage|---------------"""

C_p_25 = Perf().C_p_air(T_t_25)
gamma_air_25 = Perf().gamma(C_p_25, R)

A_25 = Gen().A(r_tip_HPC[1], r_hub_HPC[1])
mass_flow_bypass = Gen().mass_flow_bypass(mass_flow, BPR)
MFP_25 = mass_flow_bypass * np.sqrt(T_t_25) / (P_t_25 * A_25)
M_25 = Gen().M_solve(gamma_air_25, MFP_25)
T_25 = T_t_25 / Ex().stagnation_TR(gamma_air, M_2)

a_25 = Ex().a(gamma_air_25, R, T_25)
c_air_HPC = M_ax_HPC * a_25

U_HPC = DT2().U(M_rel_HPC, a_25, c_air_HPC)
omega_HPC = U_HPC / r_tip_HPC[1]
RPM_HPC = Gen().RPM(omega_HPC)

U_mid_HPC = []
for i in range(1, len(r_tip_HPC)):
    U_mid_HPC.append(DT2().U_mid(r_tip_HPC[i], r_hub_HPC[i], omega_HPC[i]))

sum_U_squared_HPC = sum(u ** 2 for u in U_mid_HPC)

Delta_H_HPC = DT2().delta_H(psi_HPC, sum_U_squared_HPC)
Delta_T_HPC = Delta_H_HPC / C_p_25

T_t_3 = T_t_25 + Delta_T_HPC
HPC_PR = Ex().stagnation_PR_fcn_TR(T_t_3, T_t_25, eta_HPC, gamma_air_25)

cooling = 1 - 0.2
P_t_3 = P_t_25 * HPC_PR * cooling

# Tillägg kylning, ca 20% Ta enbart ut kylning efter HPC, sida 104-105 Mattingly.

"""---------------|Combustor|---------------"""

gamma_air_3 = gamma_air
C_p_3 = C_p_air

'''NYTT GAMMA OCH CP FÖR GAS OCH AIR'''

# m_c = 0.05 + 0.00015 * (T_6 - 1200)  # eqn (19) 8Performance

# T_3 = T_t_3 / Ex().stagnation_TR(gamma_air, M_3)
# C_p_mix = Perf().C_p_prod(T_3)
# gamma_mix_5 = Perf().gamma(C_p_mix_5, R)

'''ÄNDRA DETTA'''
gamma_mix = (gamma_air_3+gamma_gas)/2
C_p_mix = C_p_air

Delta_H_combustion = C_p_mix  # hur mkt bränsle

P_t_4 = P_t_3 * combustor_PR  # KOLLA SÅ RÄTT RATIO, tror rätt att trycket minskar?

# massflow_mixed = massflow_ny + massflow_fuel

"""---------------|HPT - 1 stage|---------------"""

gamma_mix_2 = gamma_mix
C_p_mix_4 = C_p_air
#C_p_mix_6 = Perf().C_p_prod(T_6)

omega_HPT = omega_HPC

U_HPT = r_tip_HPT * omega_LPC

sum_U_squared_HPT = U_HPT ** 2

Delta_H_HPT = DT2().delta_H(psi_HPT, sum_U_squared_HPT)
Delta_T_HPT = Delta_H_HPT / C_p_mix_4

T_t_4 = T_turbine
T_t_45 = T_t_4 + Delta_T_HPT

HPT_PR = Ex().stagnation_PR_fcn_TR(T_t_45, T_t_4, eta_HPT, gamma_mix_2)
P_t_45 = P_t_4 * HPT_PR

"""---------------|LPT - 3 stage|---------------"""

# C_p_mix_7 = Perf().C_p_prod(T_7)
# gamma_mix_7 = Perf().gamma(C_p_mix_7, R)

gamma_mix_3 = gamma_mix
R_mix = R # detta är fel
C_p_mix_45 = C_p_air

omega_LPT = omega_fan  # samma axel

U_mid_LPT = []
for i in range(1, len(r_tip_LPT)):
    U_mid_LPT.append(r_tip_LPT[i] * omega_LPT)

sum_U_squared_LPT = sum(u ** 2 for u in U_mid_LPT)

Delta_H_LPT = DT2().delta_H(psi_LPT, sum_U_squared_LPT)
Delta_T_LPT = Delta_H_LPT / C_p_mix_45

T_t_5 = T_t_45 + Delta_T_LPT

LPT_PR = Ex().stagnation_PR_fcn_TR(T_t_5, T_t_45, eta_LPT, gamma_mix_3)
P_t_5 = P_t_4 * LPT_PR

r_ratio_BPR = r_tip_fan[1] / r_tip_LPT[3]
BPR = ((np.log(1.9176 - (r_ratio_BPR * 1.25))) / (-0.2503)) - 0.6410

# print(f'BPR = {BPR}\nMassflow = {mass_flow}')

"""---------------|Lobed mixer|---------------"""
FAR = 0.025 #Detta måste definieras bättre någonstans, förmodligen i combustor-------------------------------------
mass_flow_bypass = mass_flow * BPR/(1+BPR)
mass_flow_core =  mass_flow - mass_flow_bypass
mass_flow_hot = mass_flow_core * (1+FAR)
A_16 = 0.9367 # [m^2] manuellt beräknat från bild
A_core = 0.686 # [m^2] manuellt beräknat från bild, DETTA ÄR A_6 egentligen Ska också kunna variera med att spiken flyttas ut och in
n_bypass = 1 # Term som beskriver friktionsförluster
# A_21 = np.pi*r_tip_fan(2)^2*(1-r_tip_LPC(1)/r_tip_fan(2));
# MFP_21 = massflow_bypass*sqrt(T_3)/(P_3*A_21);
MFP_16 = mass_flow_bypass*np.sqrt(T_t_21)/(P_t_21*n_bypass*A_16) # T_3 och P_3 används ty antas isentropiskt.
M_16 = fsolve(Gen().MFP(M, gamma_air, R) - MFP_16, 0)
MFP_8 = mass_flow_hot*np.sqrt(T_t_5)/(P_t_5*A_core) # Detta borde nog ej vara A_core--------------------------
M_8 = fsolve(Gen().MFP(M, gamma_mix_3, R_mix) - MFP_16, 0)
# mixern är lobbad
T_t6A = (C_p_21 * T_t_21 * mass_flow_bypass + C_p_mix_45 * T_t_5 * mass_flow_hot) / (mass_flow_bypass + mass_flow_hot)




# parameter_list = [
#
# ]
#
# name_list = [
#
# ]
# values = [(name_list[item], parameter_list[item]) for item in range(len(parameter_list))]
# for i in range(len(values)):
#     print(values[i])
