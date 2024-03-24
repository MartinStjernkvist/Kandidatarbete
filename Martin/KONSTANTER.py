"""
import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
from PLOTTING_Funktioner import Plot_functions
from FORMLER_Exempel import Ex
from FORMLER_8Performance import Perf
from KONSTANTER import *
"""
import numpy as np

# Molar_air = 0.0289

# r_tip_fan = 0.9
# r_hub_fan = r_tip_fan - 0.62
# r_splitter = r_tip_fan - 0.316
# r_hub_splitter = r_hub_fan

P_effektivitet = 0.99  # Pt2 Verklig / Pt2 teoretisk

"""
Introduktionskompendium
"""
gamma_air = 1.4
gamma_t = 1.3  # placeholder
gamma_gas = 1.333  # sida 16

g_0 = 9.81  # gravitation
rho_SL = 1.225  # kg / m^3 air sealevel

T_std = 288.15  # K
P_std = 101_325  # Pa
R = 287.05  # J/kg*K

g_c = 1
A_fan = 2.296  # m^2

"""
Mattingly'
"""
eta_fan = 0.89  # sida 107, tabell 4.4

"""
DesignTask2
"""

EIS = 2012

"""---------------|Fan|---------------"""
M_ax_fan = 0.603
hade_angle_alpha = 15
# r_tip_splitter = 0.98 * r_tip_fan

# stage loading
psi_fan = 0.6  # sida 10: mellan 0.43 och 0.65
nu_fan = 44.29 / (98.94 + np.exp((0.01850 * EIS) - 33.31))  # hub-tip ratio

"""---------------|LPC|---------------"""

psi_LPC = -8.968 + 0.004877 * EIS  # sida 11

"""---------------|HPC|---------------"""

M_rel_HPC = 1.3
M_ax_HPC = 0.482  # sida 13
r_tip_HPC_medel = 0.366  # 0.4575 approximerat [m]
psi_HPC = -5.716 + 0.00323 * EIS  # sida 13

"""
handledning
"""
M_cruise = 1.7
h_cruise = 18288

M_rel_fan = 1.5  # hitta nytt värde
kappa = 0.8

T_turbine = 1700  # trodde det var 1650 K
T_out_compressor = 1000  # turbofan

F_symphony = 160_000  # N

eta_combustor = 0.86
combustor_PR = 0.96  # 4 procent

"""
8Performance
"""
C_p_air = 1_005.7  # J / (kg * K)
C_p_gas = 1_156.9

# s_G = ? # ground roll distance
# s_B = ? # breaking distance
# s_TO = ? # takeoff distance

"""---------------|FÖR ATT BERÄKNA C_P|---------------"""
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
https://publications.lib.chalmers.se/records/fulltext/186999/186999.pdf
sida 24
"""

T_HPC_exit = 950

# eta_fan = 0.935
eta_IPC = 0.922
# eta_HPC = 0.925
eta_HPT = 0.907
eta_IPT = 0.914
eta_LPT = 0.933

"""
https://core.ac.uk/reader/70614400
"""

epsilon_IC = 1  # placeholder, effectiveness
pi_IC = 0.96

"""
approximativa värden
"""

t_R = 3  # runtime on takeoff
s_TO = 1500  # takeoff distance in ft

eta_C = 0.95
eta_burner = 0.95
eta_turbine = 0.95
eta_mech = 0.95

# placeholder values
pi_burner = 1
pi_C = 1
pi_reg = 1
epsilon_reg = 1

Q_r = 1

"""
Värden från Katarina
"""

"""---------------|Fan|---------------"""

# Entry, approximate [m]
r_tip_fan = [0] * 3
r_hub_fan = [0] * 3

r_tip_fan[0] = 1000000  # används inte i beräkningar
r_hub_fan[0] = 1000000  # används inte i beräkningar

r_tip_fan[1] = 0.9525
r_hub_fan[1] = 0.309

# Exit, approximate [m]
r_tip_fan[2] = 0.9525
r_hub_fan[2] = 0.516



"""---------------|LPC - 3 steg|---------------"""

# Three Stage Low Pressure Compressor (LPC) Input Values
eta_LPC = 0.901  # Hitta värde för trestegs kompressor

# Stage 1, approximate [m]
r_tip_LPC = [0] * 4
r_hub_LPC = [0] * 4

r_tip_LPC[1] = 0.716
r_hub_LPC[1] = 0.545

# Stage 2, approximate [m]
r_tip_LPC[2] = 0.716
r_hub_LPC[2] = 0.545

# Stage 3, approximate [m]
r_tip_LPC[3] = 0.716
r_hub_LPC[3] = 0.545

# För att jämföra med uppmätta värden
nu_LPC_entry = 0.630  # hub-tip ratio steg 1 för LPC
nu_LPC_exit = 0.819  # hub-tip ratio steg 3 för LPC

A_21 = np.pi * r_tip_LPC[1]**2 * (1 - r_hub_LPC[1] / r_tip_LPC[1])

"""---------------|MFP|---------------"""

M_2 = 0.6  # 0.6-0.8
A_2 = 2.296  # [m2] estimerat från bild

"""---------------|HPC|---------------"""

eta_HPC = 0.941

# Stage 1, approximate [m]
r_tip_HPC = [0] * 7
r_hub_HPC = [0] * 7

r_tip_HPC[0] = 1000000  # används inte i beräkningar
r_hub_HPC[0] = 1000000  # används inte i beräkningar

r_tip_HPC[1] = 0.444
r_hub_HPC[1] = 0.273

# Stage 2, approximate [m]
r_tip_HPC[2] = 0.444
r_hub_HPC[2] = 0.334

# Stage 3, approximate [m]
r_tip_HPC[3] = 0.444
r_hub_HPC[3] = 0.364

# Stage 4, approximate [m]
r_tip_HPC[4] = 0.444
r_hub_HPC[4] = 0.371

# Stage 5, approximate [m]
r_tip_HPC[5] = 0.444
r_hub_HPC[5] = 0.400

# Stage 6, approximate [m]
r_tip_HPC[6] = 0.444
r_hub_HPC[6] = 0.436

"""---------------|Combustor|---------------"""
# Combustion Input Values


"""---------------|HPT - 3 steg|---------------"""
# Single Stage High Pressure Turbine (HPT) Input Values
r_tip_HPT = 0.509
# r_hub_HPT = 0.444 OBS ÄNDRA DESSA SEN
psi_HPT = 3.000
# eta_HPT = 0.89 # OBS ÄNDRA DESSA SEN

"""---------------|LPT - 3 steg|---------------"""
# Three Stage Low-Pressure Turbine (LPT) Input Values
psi_LPT = -39.26
# eta_LPT = 0.90 # OBS ÄNDRA DESSA SEN

# Stage 1, approximate [m]
r_tip_LPT = [0] * 4

r_tip_LPT[0] = 1000000  # används inte i beräkningar

r_tip_LPT[1] = 0.629

# Stage 2, approximate [m]
r_tip_LPT[2] = 0.768

# Stage 3, approximate [m]
r_tip_LPT[3] = 0.814

"""---------------|Lobed Mixer|---------------"""
A_8 = 1.18  # m^2 beräknat från bild
eta_nozzle = 0.98 # Beräknat mha integration????


"""---------------|Bypass and massflow|---------------"""

r_ratio_BPR = r_tip_fan[1]/r_tip_LPT[3]
BPR = ((np.log(1.9176-(r_ratio_BPR*1.25)))/(-0.2503))-0.6410