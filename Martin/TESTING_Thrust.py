import numpy as np
import matplotlib.pyplot as plt
from Thrust import F_V2
from FORMLER_AppD_V2 import AppD
from FORMLER_Exempel import Ex

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

A_0_exempel = 2.1081269749763853  # from M_0 = 0.78, mass flow = 185 kg/s, h = 10668

# temporärt
gamma_g = 1.333  # Detta är temporärt, måste fixas så den är korrekt.
cp_a = gamma_a * R / (gamma_a - 1)  # Obs! Detta gäller för ideala gaser
# cp_f = 0 # specifik värmekapacitet för bränslet
cp_g = gamma_g * R / (gamma_g - 1)  # Obs! Detta gäller för ideala gaser

# FAR = 0.025  # fuel air ratio Denna finns nu med i våran funktion, behövs nog ej längre
cp_a_sls = 10_035  # [J/kgK] värmekapacitet för luft sea level standard
# cp_f = 1  # specifik värmekapacitet för bränslet

plt.figure(figsize=(10, 10))

parameter_range_M_0 = np.linspace(0, 2, 100)
parameter_range_h = np.linspace(0, 11000, 5)

h_exempel = 10668
for h in parameter_range_h:

    M_ref = 0.78
    T_t_4_max = 1830
    alpha_list = []
    F_h_list = []
    for M in parameter_range_M_0:
        F_SL = F_V2(overall_PR, fan_PR, bypass_PR, HPC_PR, turbine_INT, fan_IPC_HPC_poly, HPT_poly, LPT_poly,
                    combustor_PL, A_0_exempel, M, h=0)
        F_h = F_V2(overall_PR, fan_PR, bypass_PR, HPC_PR, turbine_INT, fan_IPC_HPC_poly, HPT_poly, LPT_poly,
                   combustor_PL, A_0_exempel, M, h)
        F_h_list.append(F_h)
        alpha = F_h / F_SL
        theta_0 = AppD().theta_0(Ex().temperature(h), gamma_a, M)
        # TR = T_t_4_max * theta_0 / turbine_INT  # Equation (D.5) Appendix D
        # print(f'throttle ratio: {TR}')
        alpha_list.append(alpha)
    plt.plot(parameter_range_M_0, alpha_list, label=f'{str(h)} km')

plt.title('Installed Full Throttle Thrust Lapse', weight='bold')
plt.xlabel('Mach number')
plt.ylabel('Thrust Lapse (\u03B1)')
plt.legend()
plt.show()
