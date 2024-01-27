import numpy as np

g0 = 9.18
rho_SL = 1


class Master_eqn():
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
        self.P_s = dh_dt + (V / g0) * dV_dt



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
    """

    def __init__(self, T_SL, W_TO, beta, alpha, q, S, C_D0, C_DR, V, K1, K2,
                 dh_dt=0, dV_dt=0, n=1
                 ):
        Master_eqn.__init__(self, T_SL, W_TO, beta, alpha, q, S, n, C_D0, C_DR, V, K1, K2, dh_dt, dV_dt)

    def wing_loading_min(self):
        return (self.q / self.beta) * np.sqrt((self.C_D0 + self.C_DR) / self.K1)

    def thrust_to_weight_min(self):
        return ((self.n * self.beta) / self.alpha) * (
                2 * np.sqrt((self.C_D0 + self.C_DR) * self.K1) + self.K2)


class Case2(Master_eqn, Case1):
    """
    Constant Speed Climb (P_s = dh/dt)
    Given: dV/dt=0, n=1, h, dh/dt>0, V, q

    Extra:
    Known: dV_dt=0, n=1

    same wing_loading_min, thrust_to_weight_min as Case1
    """

    def __init__(self, T_SL, W_TO, beta, alpha, q, S, C_D0, C_DR, V, K1, K2, dh_dt,
                 dV_dt=0, n=1
                 ):
        Master_eqn.__init__(self, T_SL, W_TO, beta, alpha, q, S, n, C_D0, C_DR, V, K1, K2, dh_dt, dV_dt)

    def thrust_to_weight_min(self):
        return (self.beta / self.alpha) * (
                2 * np.sqrt((self.C_D0 + self.C_DR) * self.K1) + self.K2 + (self.P_s / self.V))


class Case3(Master_eqn, Case1):
    """
    Constant Altitude / Speed Turn (P_s = 0)
    Given: dh/dt=0, dV/dt=0, h, V, q, n>1

    Extra: R
    Known: dh_dt=0, dV_dt=0

    same wing_loading_min, thrust_to_weight_min as Case1
    """

    def __init__(self, T_SL, W_TO, beta, alpha, q, S, n, C_D0, C_DR, V, K1, K2,
                 R,
                 dh_dt=0, dV_dt=0
                 ):
        Master_eqn.__init__(self, T_SL, W_TO, beta, alpha, q, S, n, C_D0, C_DR, V, K1, K2, dh_dt, dV_dt)
        self.R_C = R
        self.n = np.sqrt(1 + (self.V ** 2 / (g0 ** 2 * self.R_C)) ** 2)


class Case4(Master_eqn):
    """
    Horizontal Acceleration [P_s = (V/g0)(dV/dt)]
    Given: dh/dt=0, n=1, h, V_initial, V_final, delta_t_allowable

    Known: dh_dt=0, n=1
    """

    def __init__(self, T_SL, W_TO, beta, alpha, q, S, C_D0, C_DR, V, K1, K2,
                 V_inital, V_final, delta_t_allowable, dV_dt,
                 dh_dt=0, n=1
                 ):
        Master_eqn.__init__(self, T_SL, W_TO, beta, alpha, q, S, n, C_D0, C_DR, V, K1, K2, dh_dt, dV_dt)
        self.V_initial = V_inital
        self.V_final = V_final
        self.delta_t_allowable = delta_t_allowable
    # ????????????????
    # ????????????????
    # ????????????????


class Case5(Master_eqn):
    """
    Takeoff Ground Roll (s_G), when T_SL >> (D + R)
    Given: dh/dt=0, s_G, rho, C_Lmax, V_TO = k_TO V_STALL

    Extra (parameters not in Master_eqn): C_Lmax, rho, s_G, k_TO, V_STALL
    Known: dh_dt=0
    """

    def __init__(self, T_SL, W_TO, beta, alpha, q, S, n, C_D0, C_DR, V, K1, K2, dV_dt,
                 C_Lmax, rho, s_G, k_TO, V_STALL,
                 dh_dt=0
                 ):
        Master_eqn.__init__(self, T_SL, W_TO, beta, alpha, q, S, n, C_D0, C_DR, V, K1, K2, dh_dt, dV_dt)
        self.C_Lmax = C_Lmax
        self.s_G = s_G
        self.rho = rho
        self.k_TO = k_TO
        self.V_STALL = V_STALL
        self.V_TO = k_TO * V_STALL

    def thrust_to_weight(self):
        """
        Equation (2.22)
        """
        return ((self.beta ** 2 / self.alpha) *
                (self.k_TO ** 2 / (self.s_G * self.rho * g0 * self.C_Lmax)) *
                (self.W_TO / self.S))


class Case6(Case5):
    """
    Takeoff Ground Roll (s_G)
    Given: dh/dt=0, rho, D = q C_D S, C_Lmax, V_TO = k_TO V_STALL
    and R=q C_DR S + mu_TO (beta W_TO - q C_L S)

    Extra: C_D, C_L, C_Lmax, rho, s_G, k_TO, V_STALL, mu_TO
    Known: dh_dt=0

    thrust_to_weight same as Case5
    """

    def __init__(self, T_SL, W_TO, beta, alpha, q, S, n, C_D0, C_DR, V, K1, K2, dV_dt,
                 C_D, C_L, C_Lmax, rho, s_G, k_TO, V_STALL, mu_TO,  # s_G Not known
                 dh_dt=0
                 ):
        Case5.__init__(self, T_SL, W_TO, beta, alpha, q, S, n, C_D0, C_DR, V, K1, K2, dh_dt, dV_dt,
                       C_Lmax, rho, s_G, k_TO, V_STALL)
        self.C_D = C_D
        self.C_L = C_L
        self.D = q * C_D
        self.rho = rho
        self.mu_TO = mu_TO
        self.k_TO = k_TO
        self.V_STALL = V_STALL
        self.V_TO = k_TO * V_STALL
        self.R = q * C_DR * S + mu_TO * (beta * W_TO - q * C_L * S)

    def evaluate_s_G(self):
        """
        Look at "note" page 30-31
        """
        pass


class Case7(Master_eqn):
    """
    Breaking Roll (s_B)
    Given: alpha =< 0 (reverse thrust), dh/dt=0, rho, V_TD = k_TD V_STALL, D = q C_D S
    and R=q C_DR S + mu_TO (beta W_TO - q C_L S)

    Extra: C_D, rho, s_G, k_TO, V_STALL, mu_TO
    Known: dh_dt=0
    """

    def __init__(self, T_SL, W_TO, beta, alpha, q, S, n, C_D0, C_DR, V, K1, K2, dV_dt,
                 C_Lmax, C_D, C_L, rho, k_TD, V_STALL, mu_TO,
                 dh_dt=0
                 ):
        Master_eqn.__init__(self, T_SL, W_TO, beta, alpha, q, S, n, C_D0, C_DR, V, K1, K2, dh_dt, dV_dt)
        self.C_Lmax = C_Lmax
        self.C_D = C_D
        self.C_L = C_L
        self.D = q * C_D
        self.mu_TO = mu_TO
        self.R = q * C_DR * S + mu_TO * (beta * W_TO - q * C_L * S)
        self.rho = rho
        self.k_TD = k_TD
        self.V_TD = k_TD * V_STALL

    def evaluate_s_B(self):
        pass

    def thrust_to_weight(self):
        """
        Equation (2.34)
        """
        return ((self.beta ** 2 / (- self.alpha)) *
                (self.k_TD ** 2 / (self.s_B * self.rho * g0 * self.C_Lmax)) *
                (self.W_TO / self.S))
    # ????????????????
    # ????????????????
    # ????????????????


class Case8(Master_eqn):
    """
    Service Ceiling (P_s = dh/dt)
    Given: dV/dt=0, n=1, h, dh/dt>0, C_L

    Extra: C_L
    Known: dV_dt=0, n=1
    """

    def __init__(self, T_SL, W_TO, beta, alpha, q, S, C_D0, C_DR, V, K1, K2, dh_dt,
                 C_L,
                 dV_dt=0, n=1
                 ):
        Master_eqn.__init__(self, T_SL, W_TO, beta, alpha, q, S, n, C_D0, C_DR, V, K1, K2, dh_dt, dV_dt)
        self.C_L = C_L

    def thrust_to_weight(self):
        return ((self.beta / self.alpha) *
                (self.K1 * self.C_L + self.K2 + ((self.C_D0 + self.C_DR) / self.C_L) +
                 (1 / self.V) * self.dh_dt))


class Case9(Master_eqn):
    """
    Takeoff Climb Angle
    Given: theta, n=1, dV/dt=0, C_DR, C_Lmax, k_TO, h, sigma

    Extra: theta, C_Lmax, k_TO, sigma
    Known: n=1, dV_dt=0
    """

    def __init__(self, T_SL, W_TO, beta, alpha, q, S, C_D0, C_DR, V, K1, K2, dh_dt,
                 theta, C_Lmax, k_TO, sigma,
                 n=1, dV_dt=0
                 ):
        Master_eqn.__init__(self, T_SL, W_TO, beta, alpha, q, S, n, C_D0, C_DR, V, K1, K2, dh_dt, dV_dt)
        self.theta = theta
        self.C_Lmax = C_Lmax
        self.k_TO = k_TO
        self.sigma = sigma

    def thrust_to_weight(self):
        """
        Equation (2.43)
        """
        return ((self.beta / self.alpha) *
                ((self.K1 * self.C_Lmax / self.k_TO ** 2) + self.K2 +
                 ((self.C_D0 + self.C_DR) / (self.C_Lmax / self.k_TO ** 2)) + np.sin(self.theta)))

    def evaluate_V(self):
        """
        Equation (2.44)
        """
        return np.sqrt(((2 * self.beta * self.k_TO ** 2) / (self.sigma * rho_SL * self.C_Lmax)) * (self.W_TO / self.S))



# Example

