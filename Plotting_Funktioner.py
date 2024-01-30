import numpy as np
import matplotlib.pyplot as plt

from Formler_Kap2 import MasterEqn, Case1, Case2, Case3, Case4, Case5, Case6, Case7, Case8, Case9
from Formler_Kap2 import Mach_vs_ThrustLapse


def plot_case1_instances(parameter_range, T_SL, alpha, q, S, V, dh_dt=0, dV_dt=0, n=1):
    ratios_T_W = []
    ratios_W_S = []

    for param_value in parameter_range:
        # Adjust the values of C_D0, C_DR, and K1 as needed
        case_instance = Case1(T_SL, param_value, 0.98, alpha, q, S, param_value, param_value, V, param_value,
                              param_value, dh_dt, dV_dt, n)
        ratios_T_W.append(case_instance.master_thrust_to_weight())
        ratios_W_S.append(param_value / S)

    plt.plot(ratios_W_S, ratios_T_W, label=f'Parameter Value = {param_value}')


def plot_Mach_vs_ThrustLapse(parameter_range, h, TR):
    ThrustLapse_max = []
    ThrustLapse_military = []

    for param_value in parameter_range:
        instance = Mach_vs_ThrustLapse(h=h, M_0=param_value, TR=TR)
        maximum_alpha = instance.maximum_alpha()
        military_alpha = instance.military_alpha()
        print(f'max alpha: {maximum_alpha}, mil alpha: {military_alpha}')
        ThrustLapse_max.append(maximum_alpha)
        ThrustLapse_military.append(military_alpha)

    plt.plot(parameter_range, ThrustLapse_max, label=f'Max, {TR}, h: {h}')
    plt.plot(parameter_range, ThrustLapse_military, label=f'Mil, {TR}, h: {h}')
