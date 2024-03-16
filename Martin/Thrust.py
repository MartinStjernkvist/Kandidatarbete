import numpy as np
#from Kandidatarbete.Martin.FORMLER_Exempel import Ex
from FORM_Exempel import Ex
from KONSTANTER import *
#denna används för att jag ska kunna köra filerna lokalt/ Filip


# exempelupg
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
gamma_g = 1.333  # Detta är temporärt, måste fixas så den är korrekt.
cp_a = gamma_air * R / (gamma_air - 1)  # Obs! Detta gäller för ideala gaser
# cp_f = 0 # specifik värmekapacitet för bränslet
cp_g = gamma_g * R / (gamma_g - 1)  # Obs! Detta gäller för ideala gaser

FAR = 0.025  # fuel air ratio DENNA HAR VI EN BÄTTRE APROXIMATIV FORMEL FÖR NU SOM FINNS I F FUNKTIONEN
cp_a_sls = 10_035  # [J/kgK] värmekapacitet för luft sea level standard
# cp_f = 1  # specifik värmekapacitet för bränslet

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
    a_0 = Ex().a(gamma_air, R, T_0)
    rho_0 = Ex().rho(p_0, R, T_0)
    c_0 = M_0 * a_0
    A_0 = tot_AMF / (c_0 * rho_0)

    p_02 = Ex().p_0(p_0, gamma_air, M_0)
    T_02 = Ex().T_0(T_0, gamma_air, M_0)

    p_021 = p_02 * fan_PR
    T_021 = T_02 * (p_021 / p_02) ** ((1 / eta_p_fan) * ((gamma_air - 1) / gamma_air))

    # intermediate pressure compressor
    IPC_PR = overall_PR / (fan_PR * HPC_PR)

    p_026 = p_021 * IPC_PR
    T_026 = T_021 * (p_026 / p_021) ** ((1 / eta_p_IPC) * ((gamma_air - 1) / gamma_air))

    p_03 = p_026 * HPC_PR
    T_03 = T_026 * (p_03 / p_026) ** ((1 / eta_p_HPC) * ((gamma_air - 1) / gamma_air))

    p_04 = p_03 * (1 - combustor_PL)
    T_04 = turbine_INT

    FAR = 0.5325 - 2207.5/(8200-T_03) # Detta är från ett flertal linjära approximeringar för FAR för turbine_int = 1650K. Är säkert en del fel men är bättre än fast 0.025

    mdot = tot_AMF
    mdot_bypass = mdot / (1 + (1 / bypass_PR))
    mdot_core = mdot - mdot_bypass
    mdot_fuel = mdot_core * FAR

    T_045 = T_04 - (cp_a * (T_03 - T_026)) / ((1 + FAR) * cp_g)
    p_045 = p_04 * (T_045 / T_04) ** (gamma_g / ((gamma_g - 1) * eta_p_HPT))

    Cp_a = cp_a*mdot_core
    Cp_g = cp_g*mdot_core

    T_05 = T_045 - (Cp_a / (Cp_g * (1 + FAR))) * ((1 + bypass_PR) * (T_021 - T_02) + (T_026 - T_021))
    p_05 = p_045 * (T_05 / T_045) ** (gamma_g / ((gamma_g - 1) * eta_p_LPT))

    p_018_p_18_PR = ((1 + gamma_air) / 2) ** (gamma_air / (gamma_air - 1))
    p_08_p8_PR = ((1 + gamma_g) / 2) ** (gamma_g / (gamma_g - 1))

    if (p_021 / p_0) > p_018_p_18_PR:
        bypass_choked = True
        M_18 = 1
        p_018 = p_021  # approx
        p_18 = p_018 / p_018_p_18_PR
        T_18 = Ex().T(T_021, gamma_air, M_18)
        a_18 = Ex().a(gamma_air, R, T_18)
        c_18 = M_18 * a_18
        rho_18 = Ex().rho(p_18, R, T_18)
        A_18 = mdot_bypass / (c_18 * rho_18)
        bypass_Thrust_choked = A_18 * (p_18 - p_0)

    else:
        bypass_choked = False
        p_18 = p_0
        M_18 = np.sqrt(((p_021 / p_0) ** ((gamma_air - 1) / gamma_air) - 1) * (2 / (gamma_air - 1)))  # Introkompendie sida 16
        T_18 = Ex().T(T_021, gamma_air, M_18)
        a_18 = Ex().a(gamma_air, R, T_18)
        c_18 = M_18 * a_18
        rho_18 = Ex().rho(p_18, R, T_18)
        A_18 = mdot_bypass / (c_18 * rho_18)
        bypass_Thrust_choked = 0

    if (p_05 / p_0) > p_08_p8_PR:
        core_choked = True
        M_8 = 1
        p_08 = p_05  # approx
        T_08 = T_05  # approx
        p_8 = Ex().p(p_08, gamma_g, M_8)
        T_8 = Ex().T(T_08, gamma_g, M_8)
        a_8 = Ex().a(gamma_g, R, T_8)
        c_8 = M_8 * a_8
        rho_8 = Ex().rho(p_8, R, T_8)
        A_8 = (mdot_core + mdot_fuel) / (c_8 * rho_8)
        core_Thrust_choked = A_8 * (p_8 - p_0)

    else:
        core_choked = False
        p_8 = p_0

        M_8_squared = ((p_05 / p_0) ** ((gamma_g - 1) / gamma_g) - 1) * (2 / (gamma_g - 1))  # Introkompendie sida 16
        print(f'\nM_8_squared: {M_8_squared}')
        M_8 = np.sqrt(M_8_squared)

        p_08 = Ex().p_0(p_8, gamma_g, M_8)
        T_8 = Ex().T(T_05, gamma_g, M_8)
        a_8 = Ex().a(gamma_g, R, T_8)
        c_8 = M_8 * a_8
        rho_8 = Ex().rho(p_8, R, T_8)
        A_8 = mdot_bypass / (c_8 * rho_8)
        core_Thrust_choked = 0

    F = (mdot_core + mdot_fuel) * c_8 + core_Thrust_choked + mdot_bypass * c_18 + bypass_Thrust_choked - mdot * c_0

    my_list = [p_0, T_0, a_0, rho_0, c_0, A_0, p_02, T_02, p_021, T_021, IPC_PR, p_026, T_026, p_03, T_03, p_04, T_04,
               mdot, mdot_bypass, mdot_core, mdot_fuel, T_045, p_045, T_05, p_05, p_018_p_18_PR, p_08_p8_PR, p_18, T_18,
               a_18, c_18, rho_18, A_18, p_8, T_8, a_8, c_8, rho_8, A_8]
    my_names_list = ['p_0', 'T_0', 'a_0', 'rho_0', 'c_0', 'A_0', 'p_02', 'T_02', 'p_021', 'T_021', 'IPC_PR', 'p_026',
                     'T_026', 'p_03', 'T_03', 'p_04', 'T_04', 'mdot', 'mdot_bypass', 'mdot_core', 'mdot_fuel', 'T_045',
                     'p_045', 'T_05', 'p_05', 'p_018_p_18_PR', 'p_08_p8_PR', 'p_18', 'T_18', 'a_18', 'c_18', 'rho_18',
                     'A_18', 'p_8', 'T_8', 'a_8', 'c_8', 'rho_8', 'A_8']
    values = [(my_names_list[item], my_list[item]) for item in range(len(my_list))]

    """
    print statements
    """
    print(f'\nM_0: {M_0}, h: {h}')

    for i in range(len(values)):
        # print('\n')
        print(values[i])

    print(f'\nbypass choked: {bypass_choked}')
    print(f'how choked: {bypass_Thrust_choked}')
    print(f'core choked: {core_choked}')
    print(f'how choked: {core_Thrust_choked}')
    return F


A_0_exempel = np.pi*(2.06/2)**2*(1-(2.2/9.4)**2)  # estimerad från bild


def F_V2(overall_PR, fan_PR, bypass_PR,
         HPC_PR, turbine_INT, fan_IPC_HPC_poly,
         HPT_poly, LPT_poly, combustor_PL, A_0, M_0, h):
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
    :param A:               area intake
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
    a_0 = Ex().a(gamma_air, R, T_0)
    rho_0 = Ex().rho(p_0, R, T_0)
    c_0 = M_0 * a_0
    # A_0 = tot_AMF / (c_0 * rho_0)
    A_0 = A_0

    p_02 = Ex().p_0(p_0, gamma_air, M_0)
    T_02 = Ex().T_0(T_0, gamma_air, M_0)

    p_021 = p_02 * fan_PR
    T_021 = T_02 * (p_021 / p_02) ** ((1 / eta_p_fan) * ((gamma_air - 1) / gamma_air))

    # intermediate pressure compressor
    IPC_PR = overall_PR / (fan_PR * HPC_PR)

    p_026 = p_021 * IPC_PR
    T_026 = T_021 * (p_026 / p_021) ** ((1 / eta_p_IPC) * ((gamma_air - 1) / gamma_air))

    p_03 = p_026 * HPC_PR
    T_03 = T_026 * (p_03 / p_026) ** ((1 / eta_p_HPC) * ((gamma_air - 1) / gamma_air))

    p_04 = p_03 * (1 - combustor_PL)
    T_04 = turbine_INT

    FAR = -2.77785275e-05*T_03+4.76716674e-02 #0.5325 - 2207.5/(8200-T_03) #Detta är från ett flertal linjära approximeringar för FAR för turbine_int = 1650K. Är säkert en del fel men är bättre än fast 0.025
    # mdot = tot_AMF
    M_2 = M_0 *0.28+0.2 #helt påhittat, konstig linjär aprox av M2 givet M0 för ex upgiften
    MFP = M_2 * np.sqrt(gamma_air / R) * (1 + (gamma_air - 1) / 2 * M_2 ** 2) ** ((gamma_air + 1) / (2 * (1 - gamma_air))) # 1.3 mattingly
    mdot = A_0*(p_02*MFP)/np.sqrt(T_02) # 1.3 mattingly
    mdot_bypass = mdot / (1 + (1 / bypass_PR))
    mdot_core = mdot - mdot_bypass
    mdot_fuel = mdot_core * FAR

    T_045 = T_04 - (cp_a * (T_03 - T_026)) / ((1 + FAR) * cp_g)
    p_045 = p_04 * (T_045 / T_04) ** (gamma_g / ((gamma_g - 1) * eta_p_HPT))

    Cp_a = cp_a*mdot_core
    Cp_g = cp_g*mdot_core

    T_05 = T_045 - (Cp_a / (Cp_g * (1 + FAR))) * ((1 + bypass_PR) * (T_021 - T_02) + (T_026 - T_021))
    p_05 = p_045 * (T_05 / T_045) ** (gamma_g / ((gamma_g - 1) * eta_p_LPT))

    p_018_p_18_PR = ((1 + gamma_air) / 2) ** (gamma_air / (gamma_air - 1))
    p_08_p8_PR = ((1 + gamma_g) / 2) ** (gamma_g / (gamma_g - 1))

    if (p_021 / p_0) > p_018_p_18_PR:
        bypass_choked = True
        M_18 = 1
        p_018 = p_021  # approx
        p_18 = p_018 / p_018_p_18_PR
        T_18 = Ex().T(T_021, gamma_air, M_18)
        a_18 = Ex().a(gamma_air, R, T_18)
        c_18 = M_18 * a_18
        rho_18 = Ex().rho(p_18, R, T_18)
        A_18 = mdot_bypass / (c_18 * rho_18)
        bypass_Thrust_choked = A_18 * (p_18 - p_0)

    else:
        bypass_choked = False
        p_18 = p_0
        M_18 = np.sqrt(((p_021 / p_0) ** ((gamma_air - 1) / gamma_air) - 1) * (2 / (gamma_air - 1)))  # Introkompendie sida 16
        T_18 = Ex().T(T_021, gamma_air, M_18)
        a_18 = Ex().a(gamma_air, R, T_18)
        c_18 = M_18 * a_18
        rho_18 = Ex().rho(p_18, R, T_18)
        A_18 = mdot_bypass / (c_18 * rho_18)
        bypass_Thrust_choked = 0

    if (p_05 / p_0) > p_08_p8_PR:
        core_choked = True
        M_8 = 1
        p_08 = p_05  # approx
        T_08 = T_05  # approx
        p_8 = Ex().p(p_08, gamma_g, M_8)
        T_8 = Ex().T(T_08, gamma_g, M_8)
        a_8 = Ex().a(gamma_g, R, T_8)
        c_8 = M_8 * a_8
        rho_8 = Ex().rho(p_8, R, T_8)
        A_8 = (mdot_core + mdot_fuel) / (c_8 * rho_8)
        core_Thrust_choked = A_8 * (p_8 - p_0)

    else:
        core_choked = False
        p_8 = p_0

        M_8_squared = ((p_05 / p_0) ** ((gamma_g - 1) / gamma_g) - 1) * (2 / (gamma_g - 1))  # Introkompendie sida 16
        print(f'M_8_squared: {M_8_squared}')
        M_8 = np.sqrt(M_8_squared)

        p_08 = Ex().p_0(p_8, gamma_g, M_8)
        T_8 = Ex().T(T_05, gamma_g, M_8)
        a_8 = Ex().a(gamma_g, R, T_8)
        c_8 = M_8 * a_8
        rho_8 = Ex().rho(p_8, R, T_8)
        A_8 = mdot_bypass / (c_8 * rho_8)
        core_Thrust_choked = 0

    F = (mdot_core + mdot_fuel) * c_8 + core_Thrust_choked + mdot_bypass * c_18 + bypass_Thrust_choked - mdot * c_0

    my_list = [p_0, T_0, a_0, rho_0, c_0, A_0, p_02, T_02, p_021, T_021, IPC_PR, p_026, T_026, p_03, T_03, p_04, T_04,
               mdot, mdot_bypass, mdot_core, mdot_fuel, T_045, p_045, T_05, p_05, p_018_p_18_PR, p_08_p8_PR, p_18, T_18,
               a_18, c_18, rho_18, A_18, p_8, T_8, a_8, c_8, rho_8, A_8, FAR]
    my_names_list = ['p_0', 'T_0', 'a_0', 'rho_0', 'c_0', 'A_0', 'p_02', 'T_02', 'p_021', 'T_021', 'IPC_PR', 'p_026',
                     'T_026', 'p_03', 'T_03', 'p_04', 'T_04', 'mdot', 'mdot_bypass', 'mdot_core', 'mdot_fuel', 'T_045',
                     'p_045', 'T_05', 'p_05', 'p_018_p_18_PR', 'p_08_p8_PR', 'p_18', 'T_18', 'a_18', 'c_18', 'rho_18',
                     'A_18', 'p_8', 'T_8', 'a_8', 'c_8', 'rho_8', 'A_8', 'FAR']
    values = [(my_names_list[item], my_list[item]) for item in range(len(my_list))]

    """
    print statements
    """
    print(f'\nM_0: {M_0}, h: {h}')

    for i in range(len(values)):
        # print('\n')
        print(values[i])

    print(f'\nbypass choked: {bypass_choked}')
    print(f'how choked: {bypass_Thrust_choked}')
    print(f'core choked: {core_choked}')
    print(f'how choked: {core_Thrust_choked}')
    return F


def SFC(F, tot_AMF): # detta använder det felaktiga globala FAR variabeln
    mdot = tot_AMF
    mdot_bypass = mdot / (1 + (1 / bypass_PR))
    mdot_core = mdot - mdot_bypass
    mdot_fuel = mdot_core * FAR
    return mdot_fuel / F


thrust = F(overall_PR, fan_PR, bypass_PR, HPC_PR, turbine_INT, fan_IPC_HPC_poly, HPT_poly, LPT_poly, combustor_PL,
           tot_AMF, M_0, h)
SFC = SFC(thrust, tot_AMF)
print(f'\nF (N): {thrust}')
print(f'SFC (mg/Ns): {SFC * 10 ** 6}')

M_test = 1
F_SL = F(overall_PR, fan_PR, bypass_PR, HPC_PR, turbine_INT, fan_IPC_HPC_poly, HPT_poly, LPT_poly, combustor_PL,
         tot_AMF, M_test, h=0)
print(f'F_SL: {F_SL}')

F_V2_SL = F_V2(overall_PR, fan_PR, bypass_PR, HPC_PR, turbine_INT, fan_IPC_HPC_poly, HPT_poly, LPT_poly, combustor_PL,
               A_0_exempel, M_test, h=0)
print(f'F_V2_SL: {F_V2_SL}')
