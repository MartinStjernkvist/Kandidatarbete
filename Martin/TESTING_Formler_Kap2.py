import numpy as np
import matplotlib.pyplot as plt
from Kandidatarbete.Martin.PLOTTING_Funktioner import Plot_functions
from KONSTANTER import *

gamma = gamma_a

plt.figure(figsize=(10, 8))

# Parameters for Mach_vs_ThrustLapse
parameter_range = np.linspace(0, 2, 100)

TR = [1, 1.05, 1.08]
h = 12192
symbol = ['-', '--', ':']
Plot_functions.plot_Mach_vs_ThrustLapse(parameter_range, h, gamma, TR, symbol)

TR = [1, 1.05, 1.08]
h = 0
symbol = ['-', '--', ':']
Plot_functions.plot_Mach_vs_ThrustLapse(parameter_range, h, gamma, TR, symbol)

plt.xlabel('Mach Number')
plt.ylabel('Thrust Lapse')
plt.legend()
plt.ylim((0, 1.5))

plt.show()


# Figur 2.E2
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
s_G = 1000
k_TO = 1.2

Plot_functions().plot_case5_thrust_loading_vs_wing_loading(alpha, beta, C_Lmax, rho, s_G, k_TO, wing_loading)

plt.xlabel('Wing loading')
plt.ylabel('Thrust loading')

plt.show()





