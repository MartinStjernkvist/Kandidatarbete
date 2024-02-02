import numpy as np
import matplotlib.pyplot as plt

from Plotting_Funktioner import plot_case1_instances, plot_Mach_vs_ThrustLapse
from Formler_Kap2 import MasterEqn, Case1, Case2, Case3, Case4, Case5, Case6, Case7, Case8, Case9
from Formler_Kap2 import Mach_vs_ThrustLapse
from Formler_AppD_V2 import AppendixD

# Example
# Case1_in1 = Case1(...in1 parameters...) - "in" for instance
# result_in1 = Case1_in1.thrust_to_weight() - call thrust_to_weight from Case1 onto Case1_in1
# Case1_in2 = Case1(...in2 parameters...)
# result_i2 = Case1_in2.thrust_to_weight()
# then plot result_in1 vs result_in2


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


def C_D(K1, K2, C_L, C_D0):
    return K1 * C_L ** 2 + K2 + C_D0


# Case1_in1 = Case1(T_SL=1, W_TO=1, beta=0.98, alpha=1, q=1, S=20, C_D0=1, C_DR=0, V=300, K1=1, K2=0)
# result_in1 = Case1_in1.thrust_to_weight_min()
# print(result_in1)

'''
# Parameters for Case1
T_SL = 1
alpha = 1
q = 1
S = 20
V = 300

# Create a range of parameter values (C_D0, C_DR, K1, etc.)
parameter_range = np.linspace(1, 10, 100)

plt.figure(figsize=(10, 6))
plot_case1_instances(parameter_range, T_SL, alpha, q, S, V)

plt.xlabel('W_TO/S')
plt.ylabel('T_SL/W_TO')
plt.title('Thrust-to-Weight Ratio vs. Wing Loading Ratio for Case1')
plt.legend()

plt.show()
'''

plt.figure(figsize=(10, 8))

# Parameters for Mach_vs_ThrustLapse
parameter_range = np.linspace(0, 2, 100)


TR = 1
h = 12192
symbol = '-'
plot_Mach_vs_ThrustLapse(parameter_range, h, TR, symbol)

TR = 1.05
h = 12192
symbol = '--'
plot_Mach_vs_ThrustLapse(parameter_range, h, TR, symbol)

TR = 1.08
h = 12192
symbol = ':'
plot_Mach_vs_ThrustLapse(parameter_range, h, TR, symbol)

TR = 1
h = 0
symbol = '-'
plot_Mach_vs_ThrustLapse(parameter_range, h, TR, symbol)

TR = 1.05
h = 0
symbol = '--'
plot_Mach_vs_ThrustLapse(parameter_range, h, TR, symbol)

TR = 1.08
h = 0
symbol = ':'
plot_Mach_vs_ThrustLapse(parameter_range, h, TR, symbol)


plt.xlabel('Mach Number')
plt.ylabel('Thrust Lapse')
plt.legend()
plt.ylim((0, 1.5))

plt.show()
