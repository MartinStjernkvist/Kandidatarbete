import numpy as np
import matplotlib.pyplot as plt
from Thrust import F_V2

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

# exempelupg
gamma_a = gamma
overall_PR = 45
fan_PR = 1.48
bypass_PR = 12.5
HPC_PR = 10
turbine_INT = 1650  # K
fan_IPC_HPC_poly = 0.9
HPT_poly = 0.85
LPT_poly = 0.9
combustor_PL = 0.04
tot_AMF = 185  # kg/s
A_0 = 2.1081269749763853  # from M_0 = 0.78, mass flow = 185 kg/s, h = 10668


# temporärt
gamma_g = gamma_a  # Detta är temporärt, måste fixas så den är korrekt.
cp_a = gamma_a * R / (gamma_a - 1)  # Obs! Detta gäller för ideala gaser
# cp_f = 0 # specifik värmekapacitet för bränslet
cp_g = gamma_g * R / (gamma_g - 1)  # Obs! Detta gäller för ideala gaser

FAR = 0.025  # fuel air ratio
cp_a_sls = 10_035  # [J/kgK] värmekapacitet för luft sea level standard
# cp_f = 1  # specifik värmekapacitet för bränslet
c_8 = 200
Cp_a = cp_a
Cp_g = cp_g

plt.figure(figsize=(10, 10))

parameter_range_M_0 = np.linspace(0, 2, 100)
parameter_range_h = np.linspace(0, 11000, 100)

h_ref = 10668
for h in parameter_range_h:

    M_ref = 0.78
    alpha_list = []
    F_h_list = []
    for M in parameter_range_M_0:
        F_ref = F_V2(overall_PR, fan_PR, bypass_PR, HPC_PR, turbine_INT, fan_IPC_HPC_poly, HPT_poly, LPT_poly, combustor_PL,
                 A_0, M_ref, h=h_ref)
        F_h = F_V2(overall_PR, fan_PR, bypass_PR, HPC_PR, turbine_INT, fan_IPC_HPC_poly, HPT_poly, LPT_poly, combustor_PL,
                A_0, M, h)
        F_h_list.append(F_h)
        alpha = F_h / F_ref
        alpha_list.append(alpha)
    # print(f'alpha: {alpha_list}')
    plt.plot(parameter_range_M_0, alpha_list)

plt.show()
