import numpy as np
import matplotlib.pyplot as plt

from Plotting_Funktioner import Plot_functions
from Formler_Kap2_V2 import MasterEqn, Case1, Case2, Case3, Case4, Case5, Case6, Case7, Case8, Case9
from Formler_Kap2_V2 import Mach_vs_ThrustLapse
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
g_0 = 9.81  #
rho_SL = 1.225  #
t_R = 3  # runtime on takeoff
s_TO = 500  #
T_std = 288.15  #
P_std = 101325  #
gamma = 1.4  #
R = 287  # J/kg*K
# s_G = ? # ground roll distance
# s_B = ? # breaking distance
# s_TO = ? # takeoff distance

"""
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
"""


from Plotting_Funktioner import plot_case1_thrust_loading_vs_wing_loading

plt.figure(figsize=(10, 10))
lbs_squareft_to_kg_squaremeter = 4.88242764
wing_loading = np.linspace(20, 120, 1000)

# supercruise
beta = 0.78
alpha = 1
q = 991
K1 = 0.27
K2 = 0
C_D0 = 0.028
C_DR = 0
V = 984

Plot_functions().plot_case1_thrust_loading_vs_wing_loading(alpha, beta, q, K1, K2, C_D0, C_DR, V, wing_loading)

# takeoff
C_Lmax = 2
rho = 0.002047
s_G =
k_TO
plot_case5

Plot_functions().plot_case5_thrust_loading_vs_wing_loading(alpha, beta, C_Lmax, rho, s_G, k_TO, param_value)

plt.xlabel('Wing loading')
plt.ylabel('Thrust loading')

plt.show()


