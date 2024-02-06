import numpy as np
from Formler_Kap2_V2 import MasterEqn, Case1, Case2, Case3, Case4, Case5, Case6, Case7, Case8, Case9
from Formler_Kap2_V2 import Mach_vs_ThrustLapse
from Formler_AppD_V2 import AppD
from Formler_Exempel import Ex

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
M_0 = 0.78
h = 10_668  # m

# temporärt
gamma_g = gamma_a # Detta är temporärt, måste fixas så den är korrekt.
cp_a = gamma_a*R/(gamma_a-1) # Obs! Detta gäller för ideala gaser
# cp_f = 0 # specifik värmekapacitet för bränslet
cp_g = gamma_g*R/(gamma_g-1 )# Obs! Detta gäller för ideala gaser

FAR = 0.025  # fuel air ratio
cp_a_sls = 10_035  # [J/kgK] värmekapacitet för luft sea level standard
# cp_f = 1  # specifik värmekapacitet för bränslet
c_8 = 200
Cp_a = cp_a
Cp_g = cp_g

def F(overall_PR, fan_PR, bypass_PR,
      HPC_PR, turbine_INT, fan_IPC_HPC_poly,
      HPT_poly, LPT_poly, combustor_PL, tot_AMF, M_0, h):
    """
    _PR = pressure ratio
    _R = ratio
    _T = temperature
    _poly = polytropic efficiency

    :param HPC_PR:          high pressure compressor
    :param turbine_INT:     inlet temperature
    :param HPT_poly:        high pressure turbine
    :param LPT_poly:        low pressure turbine
    :param combustor_PL:    pressure loss
    :param tot_AMF:         total air mass flow
    :param M_0:             mach number
    :param h:               height
    :return:                thrust
    """
    eta_p_fan = fan_IPC_HPC_poly
    eta_p_IPC = fan_IPC_HPC_poly
    eta_p_HPC = fan_IPC_HPC_poly
    eta_p_HPT = HPT_poly
    eta_p_LPT = LPT_poly

    _ = 1
    bypass_choked = _
    core_choked = _

    p_0 = Ex().pressure(h)
    T_0 = Ex().temperature(h)
    a_0 = Ex().a(gamma_a, R, T_0)
    rho_0 = Ex().rho(p_0, R, T_0)
    c_0 = M_0 * a_0
    A_0 = tot_AMF / (c_0 * rho_0)

    p_02 = Ex().p_0(p_0, gamma, M_0)
    T_02 = Ex().T_0(T_0, gamma, M_0)

    p_021 = p_02 * fan_PR
    T_021 = T_02 * (p_021 / p_02) ** ((1 / eta_p_fan) * ((gamma_a - 1) / gamma_a))

    # intermediate pressure compressor
    IPC_PR = overall_PR / (fan_PR * HPC_PR)

    p_026 = p_021 * IPC_PR
    T_026 = T_021 * (p_026 / p_021) ** ((1 / eta_p_IPC) * ((gamma_a - 1) / gamma_a))

    p_03 = p_026 * HPC_PR
    T_03 = T_026 * (p_03 / p_026) ** ((1 / eta_p_HPC) * ((gamma_a - 1) / gamma_a))

    p_04 = p_03 * (1 - combustor_PL)
    T_04 = turbine_INT

    mdot = tot_AMF
    mdot_bypass = mdot / (1 + (1 / bypass_PR))
    mdot_core = mdot - mdot_bypass
    mdot_fuel = mdot_core * FAR

    T_045 = T_04 - (cp_a * (T_03 - T_026)) / ((1 + FAR) * cp_g)
    p_045 = (T_045 / T_04) ** (gamma_g / ((gamma_g - 1) * eta_p_HPT))

    T_05 = T_045 - (Cp_a / (Cp_g * (1 + FAR))) * ((1 + bypass_PR) * (T_021 - T_02) + (T_026 - T_021))
    p_05 = p_045 * (T_05 / T_045) ** (gamma_g / ((gamma_g - 1) * eta_p_LPT))

    p_018_p_18_PR = ((1 + gamma_a) / 2) ** (gamma_a / (gamma_a - 1))
    p_08_p8_PR = ((1 + gamma_a) / 2) ** (gamma_a / (gamma_a - 1))

    if (p_021 / p_0) > p_018_p_18_PR and (p_05 / p_0) > p_08_p8_PR:

        choked = True
        print(f'choked: {choked}')
        # bypass choked
        M_18 = 1
        p_018 = p_021  # approx
        p_18 = p_018 / p_018_p_18_PR
        T_18 = Ex().T(T_021, gamma_a, 1)
        a_18 = Ex().a(gamma_a, R, T_18)
        c_18 = M_18 * a_18
        rho_18 = Ex().rho(p_18, R, T_18)
        A_18 = mdot_bypass / (c_18 * rho_18)

        # core choked
        M_8 = 1
        p_08 = p_05  # approx
        T_08 = T_05  # approx
        p_8 = Ex().p(p_08, gamma_g, 1)
        T_8 = Ex().T(T_08, gamma_g, 1)
        a_8 = Ex().a(gamma_g, R, T_8)
        c_8 = M_8 * a_8
        rho_8 = Ex().rho(p_8, R, T_8)
        A_8 = (mdot_core + mdot_fuel) / (c_8 * rho_8)
        F = (mdot_core + mdot_fuel) * c_8 + (p_8 - p_08) + mdot_bypass * c_18 + A_18 * (p_18 - p_0) - mdot * c_0
        SFC = mdot_fuel / F

        print()

        return F, SFC

    else:

        choked = False
        print(f'choked: {choked}')


print(F(overall_PR, fan_PR, bypass_PR,
        HPC_PR, turbine_INT, fan_IPC_HPC_poly,
        HPT_poly, LPT_poly, combustor_PL, tot_AMF, M_0, h))