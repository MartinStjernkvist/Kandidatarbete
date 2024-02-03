import numpy as np
from Formler_Kap2_V2 import MasterEqn, Case1, Case2, Case3, Case4, Case5, Case6, Case7, Case8, Case9
from Formler_Kap2_V2 import Mach_vs_ThrustLapse
from Formler_AppD_V2 import AppendixD
from Formler_Annat import Exempel

"""
Constants (keep track)
"""
g_0 = 9.81  #
rho_SL = 1.225  #
t_R = 3  # runtime on takeoff
s_TO = 1500  # takeoff distance in ft
T_std = 288.15  #
P_std = 101325  #
gamma = 1.4  #
R = 287.05  # J/kg*K

# s_G = ? # ground roll distance
# s_B = ? # breaking distance
# s_TO = ? # takeoff distance


theta_0_break = 1.2
T_5 = 10 ** 20  # limiting
stagnation_pressure_loss_combustor = 0.04
M_SL = 1
M_cruise = 1.7
BPR = 3  # around 3-4
eta_pc = 0.99

T_05 = Exempel().T_0(T_5, gamma, 1)
T_8 = Exempel().T(T_05, gamma, 1)

inlet_pressure_ratio_cruise = Exempel().SPR(gamma, M_cruise)
compressor_temperature_ratio_cruise = Exempel().STR_compressor(gamma, M_cruise, eta_pc)
print(inlet_pressure_ratio_cruise)
print(compressor_temperature_ratio_cruise)
