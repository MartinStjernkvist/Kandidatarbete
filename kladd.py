
"""
Constants (keep track)
"""
g_0 = 9.81          #
rho_SL = 1.225      #
t_R = 3             # runtime on takeoff
s_TO = 500          #
T_std = 288.15      #
P_std = 101325      #
gamma = 1.4         #
R = 287             # J/kg*K



def calculate_theta_0(T, M_0):
    """
    Equation (2.52a)
    """
    return (T / T_std) * (1 + ((gamma - 1) / 2) * M_0 ** 2)


def calculate_delta_0(P, M_0):
    """
    Equation (2.52b)
    """
    return (P / P_std) * (1 + ((gamma - 1) / 2) * M_0 ** 2) ** (gamma / (gamma - 1))


def calculate_temperature(h):
    """
    Introduktionskompendium: Equation (1.29a)
    """
    return T_std - (6.5 * 10 ** (-3) * h)


def calculate_pressure(h):
    """
    Introduktionskompendium: Equation (1.29b)
    """
    if h <= 11000:
        return P_std * ((T_std - (6.5 * 10 ** (-3) * h)) / T_std) ** (g_0 / (6.5 * 10 ** (-3) * R))
    else:
        return P_std * ((T_std - (6.5 * 10 ** (-3) * 11000)) / T_std) ** (g_0 / (6.5 * 10 ** (-3) * R))


# h = 12192
h = 11000
M_0 = 0
TR = 1

T = calculate_temperature(h)
P = calculate_pressure(h)
print(f'\ntemp: {T}, press: {P}')

print(f'\nheight: {h}, M0 = {M_0}, TR = {TR}\n')

print(f'temperature: {T}')
theta0 = calculate_theta_0(T, M_0)
print(f'theta0: {theta0}')
delta0 = calculate_delta_0(P, M_0)
print(f'delta0: {delta0}')

alpha = delta0 * (1 - 3.5 * (theta0 - TR) / theta0)
print(f'alpha: {alpha}')
