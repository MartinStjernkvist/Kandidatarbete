import numpy as np
import matplotlib.pyplot as plt

from Plotting_Funktioner import plot_case1_thrust_loading_vs_wing_loading
from Formler_Kap2_V2 import MasterEqn, Case1, Case2, Case3, Case4, Case5, Case6, Case7, Case8, Case9

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

plot_case1_thrust_loading_vs_wing_loading(wing_loading, beta, alpha, q, K1, K2, C_D0, C_DR, V)

plt.xlabel('Wing loading')
plt.ylabel('Thrust loading')

plt.show()
