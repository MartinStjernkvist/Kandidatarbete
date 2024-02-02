import numpy as np
from scipy import integrate

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


class CommonFunctionality:
    """
    Common functionality for MasterEqn and other classes
    """

    def P_s(self, V, dh_dt, dV_dt):
        """
        Weight specific excess power

        Equation (2.2b)

        :param V:       velocity
        :param dh_dt:   altitude derivative
        :param dV_dt:   acceleration
        """
        return dh_dt + (V / g_0) * dV_dt

    def temperature(self, h):
        """
        Introkompendium: Equation (1.29a)

        :param h:   altitude
        """
        T_0 = T_std - (6.5 * 10 ** (-3) * h)
        return T_0

    def pressure(self, h):
        """
        Introkompendium: Equation (1.29b)

        :param h:   altitude
        """
        if h <= 11000:
            P_0 = P_std * ((T_std - (6.5 * 10 ** (-3) * h)) / T_std) ** (g_0 / (6.5 * 10 ** (-3) * R))
            return P_0
        else:
            P_0 = P_std * ((T_std - (6.5 * 10 ** (-3) * 11000)) / T_std) ** (g_0 / (6.5 * 10 ** (-3) * R))
            return P_0

    def theta_0(self, T, M_0):
        """
        dimensionless ratio of freestream total temperature to
        sea level static temperature of the standard atmosphere

        Equation (2.52a)

        :param T:       thrust
        :param M_0:     mach number
        """
        return (T / T_std) * (1 + ((gamma - 1) / 2) * M_0 ** 2)

    def delta_0(self, P, M_0):
        """
        Equation (2.52b)
        """
        return (P / P_std) * (1 + ((gamma - 1) / 2) * M_0 ** 2) ** (gamma / (gamma - 1))

    def centripetal_force(self, V, R_C):
        """
        = n

        Equation (2.17)

        :param V:   velocity
        :param R_C: radius
        """
        n = np.sqrt(1 + (V ** 2 / (g_0 ** 2 * R_C)) ** 2)
        return n

    def dV_dt(self, V_final, V_initial, delta_t_allowable):
        return ((V_final - V_initial) / delta_t_allowable)

    def s_G(self, ksi_TO, beta, W_TO, S, rho, alpha, C_Lmax, k_TO):
        """
        ground roll distance

        Equation (2.25)

        Look at "note" page 30-31
        """
        pass

    def s_TO(self):
        # """
        # takeoff distance
        #
        # Equation (2.E1)
        # """
        # return ((self.k_TO ** 2 * self.beta ** 2) /
        #         (self.rho * self.C_Lmax * self.alpha * (self.T_SL / self.W_TO))) * (self.W_TO / self.S) + (
        #         t_R * self.k_TO) * np.sqrt((2 * self.beta) / (self.rho * self.C_Lmax)) * np.sqrt(self.W_TO / self.S)
        pass


class MasterEqn(CommonFunctionality):
    def __init__(self, T_SL, W_TO, beta, alpha, q, S, n, C_D0, C_DR, V, K1, K2, dh_dt, dV_dt):
        """
        :param T_SL:    thrust loading at sea level
        :param W_TO:    weight at takeoff
        :param beta:    instantaneous weight fraction
        :param alpha:   installed full throttle thrust lapse
        :param q:
        :param S:       surface of
        :param n:       load factor, number of g's
        :param C_D0:    drag coefficient at zero lift
        :param C_DR:    drag coefficient from external sources
        :param V:       velocity
        :param K1:
        :param K2:
        :param dh_dt:   altitude derivative
        :param dV_dt:   acceleration
        """
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
        self.P_s = self.P_s(V, dh_dt, dV_dt)

    def installed_thrust(self):
        """
        = T

        Equation (2.3)
        """
        return self.alpha * self.T_SL

    def instantaneous_weight(self):
        """
        = W

        Equation (2.4)
        """
        return self.beta * self.W_TO

    def master_thrust_to_weight(self):
        """
        = T_SL / W_TO

        Equation (2.11)
        """
        return (self.beta / self.alpha) * ((self.q * self.S) / (self.beta * self.W_TO) * (
                self.K1 * ((self.n * self.beta * self.W_TO) / (self.q * self.S)) ** 2
                + self.K2 * ((self.n * self.beta * self.W_TO) / (self.q * self.S))
                + self.C_D0 + self.C_DR) + (self.P_s / self.V))


class Case1(MasterEqn):
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
        """
        Corresponds to maximum range
        """
        return ((self.n * self.beta) / self.alpha) * (
                2 * np.sqrt((self.C_D0 + self.C_DR) * self.K1) + self.K2)


class Case2(Case1):
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
        MasterEqn.__init__(self, T_SL, W_TO, beta, alpha, q, S, n, C_D0, C_DR, V, K1, K2, dh_dt, dV_dt)

    def thrust_to_weight_min(self):
        return (self.beta / self.alpha) * (
                2 * np.sqrt((self.C_D0 + self.C_DR) * self.K1) + self.K2 + (self.P_s / self.V))


class Case3(Case1):
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
        """
        :param R:   radius of turn
        """
        MasterEqn.__init__(self, T_SL, W_TO, beta, alpha, q, S, n, C_D0, C_DR, V, K1, K2, dh_dt, dV_dt)
        self.R_C = R
        self.n = self.centripetal_force(V, self.R_C)


class Case4(MasterEqn):
    """
    Horizontal Acceleration [P_s = (V/g0)(dV/dt)]
    Given: dh/dt=0, n=1, h, V_initial, V_final, delta_t_allowable

    Extra: V_inital, V_final, delta_t_allowable
    Known: dh_dt=0, n=1

    Page 52: Mission phase 7-8: Horizontal acceleration: C_DR = K2 = 0
    """

    def __init__(self, T_SL, W_TO, beta, alpha, q, S, C_D0, C_DR, V, K1, K2, dV_dt,
                 V_initial, V_final, delta_t_allowable,
                 dh_dt=0, n=1
                 ):
        """
        :param V_inital:            initial velocity
        :param V_final:             final velocity
        :param delta_t_allowable:   allowable time to accelerate
        """
        MasterEqn.__init__(self, T_SL, W_TO, beta, alpha, q, S, n, C_D0, C_DR, V, K1, K2, dh_dt, dV_dt)
        self.V_initial = V_initial
        self.V_final = V_final
        self.delta_t_allowable = delta_t_allowable
        self.dV_dt = self.dV_dt(V_final, V_initial, delta_t_allowable)

    # def delta_t_function(self, V):
    #     return (1 / (2 * g_0)) * (V / self.P_s)
    #
    # def delta_t(self,):
    #     result, error = integrate.quad(delta_t_function, self.V_initial, self.V_final)
    #     delta_t_allowable = result
    #     return delta_t_allowable

    # ????????????????
    # ????????????????
    # ????????????????


class Case5(MasterEqn):
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
        """
        :param C_Lmax:  maximum coefficient of lift
        :param rho:
        :param s_G:     ground roll distance
        :param k_TO:
        :param V_STALL: stall velocity
        """
        MasterEqn.__init__(self, T_SL, W_TO, beta, alpha, q, S, n, C_D0, C_DR, V, K1, K2, dh_dt, dV_dt)
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
                (self.k_TO ** 2 / (self.s_G * self.rho * g_0 * self.C_Lmax)) *
                (self.W_TO / self.S))

    def wing_loading(self):
        """
        Equation (2.E2)
        """
        a = (self.k_TO ** 2 * self.beta ** 2) / (self.rho * self.C_Lmax * self.alpha * (self.T_SL / self.W_TO))
        b = (t_R * self.k_TO) * np.sqrt((2 * self.beta) / (self.rho * self.C_Lmax))
        c = s_TO
        return (-b * np.sqrt(b ** 2 * (4 * a * c) / (2 * a))) ** 2


class Case6(Case5):
    """
    Takeoff Ground Roll (s_G)
    Given: dh/dt=0, rho, D = q C_D S, C_Lmax, V_TO = k_TO V_STALL
    and R=q C_DR S + mu_TO (beta W_TO - q C_L S)

    Extra: C_D, C_L, C_Lmax, rho, s_G, k_TO, V_STALL, mu_TO
    Known: dh_dt=0

    thrust_to_weight same as Case5

    Page 50: Mission phase 1-2: Takeoff, no obstacle: s_TO = s_G + s_R
    """

    def __init__(self, T_SL, W_TO, beta, alpha, q, S, n, C_D0, C_DR, V, K1, K2, dV_dt,
                 C_D, C_L, C_Lmax, rho, s_G, k_TO, V_STALL, mu_TO,  # s_G Not known
                 dh_dt=0
                 ):
        """
        :param C_D:     drag coefficient
        :param C_L:     lift coefficient
        :param mu_TO:
        """
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
        self.ksi_TO = C_D + C_DR - mu_TO * C_L

    def thrust_to_weight(self):
        return


class Case7(MasterEqn):
    """
    Breaking Roll (s_B)
    Given: alpha =< 0 (reverse thrust), dh/dt=0, rho, V_TD = k_TD V_STALL, D = q C_D S
    and R=q C_DR S + mu_TO (beta W_TO - q C_L S)

    Extra: C_D, rho, s_G, k_TO, V_STALL, mu_TO
    Known: dh_dt=0

    Page 53: Mission phase 13-14: Landing, no reverse thrust: s_L = s_FR + s_B
    """

    def __init__(self, T_SL, W_TO, beta, alpha, q, S, n, C_D0, C_DR, V, K1, K2, dV_dt,
                 C_Lmax, C_D, C_L, rho, k_TD, V_STALL, mu_TO,
                 dh_dt=0
                 ):
        MasterEqn.__init__(self, T_SL, W_TO, beta, alpha, q, S, n, C_D0, C_DR, V, K1, K2, dh_dt, dV_dt)
        self.C_Lmax = C_Lmax
        self.C_D = C_D
        self.C_L = C_L
        self.D = q * C_D
        self.mu_TO = mu_TO
        self.R = q * C_DR * S + mu_TO * (beta * W_TO - q * C_L * S)
        self.rho = rho
        self.k_TD = k_TD
        self.V_TD = k_TD * V_STALL

    def thrust_to_weight(self):
        """
        Equation (2.34)
        Reverse thrust
        """
        return ((self.beta ** 2 / (- self.alpha)) *
                (self.k_TD ** 2 / (self.s_B * self.rho * g_0 * self.C_Lmax)) *
                (self.W_TO / self.S))
    # ????????????????
    # ????????????????
    # ????????????????


class Case8(MasterEqn):
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
        MasterEqn.__init__(self, T_SL, W_TO, beta, alpha, q, S, n, C_D0, C_DR, V, K1, K2, dh_dt, dV_dt)
        self.C_L = C_L

    def thrust_to_weight(self):
        return ((self.beta / self.alpha) *
                (self.K1 * self.C_L + self.K2 + ((self.C_D0 + self.C_DR) / self.C_L) +
                 (1 / self.V) * self.dh_dt))


class Case9(MasterEqn):
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
        """
        :param sigma:
        """
        MasterEqn.__init__(self, T_SL, W_TO, beta, alpha, q, S, n, C_D0, C_DR, V, K1, K2, dh_dt, dV_dt)
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

    def V(self):
        """
        Equation (2.44)
        """
        return np.sqrt(((2 * self.beta * self.k_TO ** 2) / (self.sigma * rho_SL * self.C_Lmax)) * (self.W_TO / self.S))


class Takeoff():
    pass


class Supercruise():
    pass


class Mach_vs_ThrustLapse(CommonFunctionality):
    """
    Functionality to create plot on page 42
    """

    def __init__(self, h, M_0, TR):
        """
        :param h:   altitude
        :param M_0: mach number
        :param TR:  throttle ratio
        """
        self.h = h
        self.T = self.temperature(h)
        self.P = self.pressure(h)
        self.M_0 = M_0
        self.TR = TR
        self.theta_0 = self.theta_0(self.T, self.M_0)
        self.delta_0 = self.delta_0(self.P, self.M_0)

    def alpha_max(self):
        """
        Equation (2.54a)
        """
        if self.theta_0 <= self.TR:
            return self.delta_0
        else:
            return self.delta_0 * (1 - ((3.5 * (self.theta_0 - self.TR)) / self.theta_0))

    def alpha_military(self):
        """
        Equation (2.54b)
        """
        if self.theta_0 <= self.TR:
            return 0.6 * self.delta_0
        else:
            return 0.6 * self.delta_0 * (1 - ((3.8 * (self.theta_0 - self.TR)) / self.theta_0))
