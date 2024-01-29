import numpy as np
import matplotlib.pyplot as plt

# from Formler_Kap2 import Master_eqn, Case1, Case2, Case3, Case4, Case5, Case6, Case7, Case8, Case9

import numpy as np

g0 = 9.81

class CommonFunctionality:
    """
    Common functionality for both Master_eqn and Case1
    """
    def calculate_P_s(self, V, dh_dt, dV_dt):
        return dh_dt + (V / g0) * dV_dt

class Master_eqn(CommonFunctionality):
    """
    Equation (2.11)
    """

    def __init__(self, T_SL, W_TO, beta, alpha, q, S, n, C_D0, C_DR, V, K1, K2, dh_dt, dV_dt):
        self.T_SL = T_SL
        self.W_TO = W_TO
        self.beta = beta
        self.alpha = alpha
        self.q = q
        self.S = S
        self.n = n
        self.C_D0 = C_D0
        self.C_DR = C_DR
        self.V = V
        self.K1 = K1
        self.K2 = K2
        self.dh_dt = dh_dt
        self.dV_dt = dV_dt
        self.P_s = self.calculate_P_s(V, dh_dt, dV_dt)

    def master_thrust_to_weight(self):
        return (self.beta / self.alpha) * ((self.q * self.S) / (self.beta * self.W_TO) * (
                self.K1 * ((self.n * self.beta * self.W_TO) / (self.q * self.S)) ** 2
                + self.K2 * ((self.n * self.beta * self.W_TO) / (self.q * self.S))
                + self.C_D0 + self.C_DR) + (self.P_s / self.V))

class Case1(Master_eqn):
    """
    Constant Altitude / Speed Cruise (P_s = 0)
    Given: dh/dt=0, dV/dt=0, n=1, h, V, q

    Extra:
    Known: dh_dt=0, dV_dt=0, n=1

    Page 44: Supercruise: C_DR = K2 = 0
    Page 51: Mission phases 6-7 and 8-9: Supersonic penetration and escape dash: C_DR = K2 = 0
    Page 54: Maximum Mach number: C_DR = K2 = 0
    """

    def __init__(self, T_SL, W_TO, beta, alpha, q, S, C_D0, C_DR, V, K1, K2,
                 dh_dt=0, dV_dt=0, n=1
                 ):
        super().__init__(T_SL, W_TO, beta, alpha, q, S, n, C_D0, C_DR, V, K1, K2, dh_dt, dV_dt)

    def wing_loading_min(self):
        return (self.q / self.beta) * np.sqrt((self.C_D0 + self.C_DR) / self.K1)

    def thrust_to_weight_min(self):
        return ((self.n * self.beta) / self.alpha) * (
                2 * np.sqrt((self.C_D0 + self.C_DR) * self.K1) + self.K2)

g0 = 9.18
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