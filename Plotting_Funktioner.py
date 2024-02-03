import numpy as np
import matplotlib.pyplot as plt

from Formler_Kap2_V2 import MasterEqn, Case1, Case2, Case3, Case4, Case5, Case6, Case7, Case8, Case9
from Formler_Kap2_V2 import Mach_vs_ThrustLapse


class Plot_functions:

    def plot_Mach_vs_ThrustLapse(parameter_range, h, TR, symbol):

        for i in range(3):
            ThrustLapse_max = []
            ThrustLapse_military = []
            for param_value in parameter_range:
                instance = Mach_vs_ThrustLapse(h=h, M_0=param_value, TR=TR[i])
                maximum_alpha = instance.alpha_max()
                military_alpha = instance.alpha_military()
                ThrustLapse_max.append(maximum_alpha)
                ThrustLapse_military.append(military_alpha)

            plt.plot(parameter_range, ThrustLapse_max, linestyle=symbol[i], linewidth=3, label=f'Max, {TR}, h: {h}')
            plt.plot(parameter_range, ThrustLapse_military, linestyle=symbol[i], linewidth=3, label=f'Mil, {TR}, h: {h}')

    def plot_case1_thrust_loading_vs_wing_loading(self, alpha, beta, q, K1, K2, C_D0, C_DR, V, wing_loading):
        thrust_loading = []

        for param_value in wing_loading:
            thrust_to_weight = Case1().thrust_to_weight(alpha, beta, q, K1, K2, C_D0, C_DR, V, param_value)
            thrust_loading.append(thrust_to_weight)

        plt.plot(wing_loading, thrust_loading)

    def plot_case5_thrust_loading_vs_wing_loading(self, alpha, beta, C_Lmax, rho, s_G, k_TO, wing_loading):
        thrust_loading = []

        for param_value in wing_loading:
            thrust_to_weight = Case5().thrust_to_weight(alpha, beta, C_Lmax, rho, s_G, k_TO, param_value)
            thrust_loading.append(thrust_to_weight)

        plt.plot(wing_loading, thrust_loading)
