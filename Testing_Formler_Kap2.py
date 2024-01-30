import numpy as np
import matplotlib.pyplot as plt

from Formler_Kap2 import MasterEqn, Case1, Case2, Case3, Case4, Case5, Case6, Case7, Case8, Case9
# Example
# Case1_in1 = Case1(...in1 parameters...) - "in" for instance
# result_in1 = Case1_in1.thrust_to_weight() - call thrust_to_weight from Case1 onto Case1_in1
# Case1_in2 = Case1(...in2 parameters...)
# result_i2 = Case1_in2.thrust_to_weight()
# then plot result_in1 vs result_in2



g0 = 9.81
rho_SL = 1.225
t_R = 3
s_TO = 500



def C_D(K1, K2, C_L, C_D0):
    return K1 * C_L **2 + K2 + C_D0

# Case1_in1 = Case1(T_SL=1, W_TO=1, beta=0.98, alpha=1, q=1, S=20, C_D0=1, C_DR=0, V=300, K1=1, K2=0)
# result_in1 = Case1_in1.thrust_to_weight_min()
# print(result_in1)

def plot_case1_instances(parameter_range, T_SL, alpha, q, S, V, dh_dt=0, dV_dt=0, n=1):
    ratios_T_W = []
    ratios_W_S = []

    for param_value in parameter_range:
        # Adjust the values of C_D0, C_DR, and K1 as needed
        case_instance = Case1(T_SL, param_value, 0.98, alpha, q, S, param_value, param_value, V, param_value, param_value, dh_dt, dV_dt, n)
        ratios_T_W.append(case_instance.master_thrust_to_weight())
        ratios_W_S.append(param_value / S)

    plt.plot(ratios_W_S, ratios_T_W, label=f'Parameter Value = {param_value}')

# Parameters for Case1
T_SL = 1
alpha = 1
q = 1
S = 20
V = 300

# Create a range of parameter values (C_D0, C_DR, K1, etc.)
parameter_range = np.linspace(1, 10, 100)

# Plotting
plt.figure(figsize=(10, 6))
plot_case1_instances(parameter_range, T_SL, alpha, q, S, V)

# Add labels and legend
plt.xlabel('W_TO/S')
plt.ylabel('T_SL/W_TO')
plt.title('Thrust-to-Weight Ratio vs. Wing Loading Ratio for Case1')
plt.legend()

# Show the plot
plt.show()