import numpy as np
from Formler_Kap2_V2 import MasterEqn, Case1, Case2, Case3, Case4, Case5, Case6, Case7, Case8, Case9
from Formler_Kap2_V2 import Mach_vs_ThrustLapse
from Formler_AppD_V2 import AppD
from Formler_Exempel import Ex

h = 10668
T = Ex().temperature(h)
print(T)
M_0 = 0.78
gamma = 1.4
print(Ex().T_0(T, gamma, M_0))