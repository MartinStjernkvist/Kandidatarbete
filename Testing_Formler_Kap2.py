import numpy as np
import matplotlib.pyplot as plt

from Formler_Kap2 import Master_eqn, Case1, Case2, Case3, Case4, Case5, Case6, Case7, Case8, Case9
# Example
# Case1_in1 = Case1(...in1 parameters...) - "in" for instance
# result_in1 = Case1_in1.thrust_to_weight() - call thrust_to_weight from Case1 onto Case1_in1
# Case1_in2 = Case1(...in2 parameters...)
# result_i2 = Case1_in2.thrust_to_weight()
# then plot result_in1 vs result_in2


import numpy as np



g0 = 9.81
rho_SL = 1.225
t_R = 3
s_TO = 500

def C_D(K1, K2, C_L, C_D0):
    return K1 * C_L **2 + K2 + C_D0

ax = plt.axes()

plt.show()

Case1_in1 = Case1(1,1,1,1,1,1,1,1,1,1,1,1,1,1)
result_in1 = Case1_in1.thrust_to_weight_min()
print(result_in1)