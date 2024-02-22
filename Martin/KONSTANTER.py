"""
import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
from PLOTTING_Funktioner import Plot_functions
from FORMLER_Exempel import Ex
from FORMLER_8Performance import Perf
from KONSTANTER import *
"""

"""in m"""
r_tip_fan = 0.9
r_hub_fan = r_tip_fan - 0.62
r_splitter = r_tip_fan - 0.316
r_hub_splitter = r_hub_fan

M_cruise = 1.7
h_cruise = 18288

gamma_a = 1.4
gamma_c = 1.33
gamma_t = 1.3 # placeholder

g_0 = 9.81  # gravitation
rho_SL = 1.225  # kg / m^3 air sealevel

T_std = 288.15  # K
P_std = 101_325  # Pa
R = 287.05  # J/kg*K

T_turbine = 1650  # K

g_c = 1
A_fan = 2.296  # m^2

F_symphony = 160_000  # N

P_effektivitet = 0.99  # Pt2 Verklig / Pt2 teoretisk

"""
värden från 8Performance
"""
Cp_air = 1_005.7  # J / (kg * K)
Cp_gas = 1_156.9

# s_G = ? # ground roll distance
# s_B = ? # breaking distance
# s_TO = ? # takeoff distance


"""
https://publications.lib.chalmers.se/records/fulltext/186999/186999.pdf
sida 24
"""

T_HPC_exit = 950

eta_fan = 0.935
eta_IPC = 0.922
eta_HPC = 0.925
eta_HPT = 0.907
eta_IPT = 0.914
eta_LPT = 0.933

"""
https://core.ac.uk/reader/70614400
"""

epsilon_IC = 1 # placeholder
pi_IC = 0.96

"""
8Performance
"""
A_0_air = 1047.63
A_1_air = - 0.39
A_2_air = 8.89 * 10 ** -4
A_3_air = - 1.64 * 10 ** -7
A_4_air = - 6.65 * 10 ** -10
A_5_air = 6.03 * 10 ** -13
A_6_air = -2.07 * 10 ** -16
A_7_air = 2.59 * 10 ** -20

A_0_prod = 309.08
A_1_prod = 9.24
A_2_prod = -1.87 * 10 ** -2
A_3_prod = 2.43 * 10 ** -5
A_4_prod = -1.85 * 10 ** -8
A_5_prod = 8.08 * 10 ** -12
A_6_prod = -1.9 * 10 ** -15
A_7_prod = 1.86 * 10 ** -19

"""
approximativa värden
"""

t_R = 3  # runtime on takeoff
s_TO = 1500  # takeoff distance in ft

eta_C = 0.95
eta_burner = 0.95
eta_turbine = 0.95
eta_mech = 0.95

pi_burner = 1
pi_C = 1
pi_reg = 1
epsilon_reg = 1

Q_r = 1
