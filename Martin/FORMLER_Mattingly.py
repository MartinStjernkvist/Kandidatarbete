import numpy as np
from scipy import integrate

# alpha, beta, q, n, K1, K2, C_D0, C_DR, V, dV_dt, dh_dt, C_D, C_L, C_Lmax, rho, s_G, s_B, k_TO, k_TD, mu_TO, V_STALL, sigma, theta, T_SL, W_TO, S, wing_loading

"""
Constants (keep track)
"""
g_0 = 9.81  #
rho_SL = 1.225  #
t_R = 3  # runtime on takeoff
s_TO = 1500  # takeoff distance in ft
T_std = 288.15  #
P_std = 101325  #
gamma = 1.4  #
R = 287.05  # J/kg*K


# s_G = ? # ground roll distance
# s_B = ? # breaking distance
# s_TO = ? # takeoff distance


def MFP(M, g_c=1, gamma=gamma):
    """
    mass flow parameter

    Equation (1.3)
    :param M:       mach number
    :param g_c:     gas constant?
    :param gamma:
    """
    return M * np.sqrt((gamma * g_c / R)) * (1 + ((gamma - 1) / 2) * M ** 2) ** ((gamma + 1) / (2 * (1 - gamma)))


class CommonFunctionality:
    """
    Common functionality for MasterEqn and other classes
    """

    def P_s(self, V, dh_dt, dV_dt):
        """
        Weight specific excess power

        Equation (2.2b)
        """
        return dh_dt + (V / g_0) * dV_dt

    def temperature(self, h):
        """
        Introkompendium: Equation (1.29a)
        """
        T_0 = T_std - (6.5 * 10 ** (-3) * h)
        return T_0

    def pressure(self, h):
        """
        Introkompendium: Equation (1.29b)
        """
        if h <= 11000:
            P_0 = P_std * ((T_std - (6.5 * 10 ** (-3) * h)) / T_std) ** (g_0 / (6.5 * 10 ** (-3) * R))
            return P_0
        else:
            P_0 = P_std * ((T_std - (6.5 * 10 ** (-3) * 11000)) / T_std) ** (g_0 / (6.5 * 10 ** (-3) * R))
            return P_0

    def theta_0(self, T, gamma, M_0):
        """
        dimensionless ratio of freestream total temperature to
        sea level static temperature of the standard atmosphere

        Equation (2.52a)
        """
        return (T / T_std) * (1 + ((gamma - 1) / 2) * M_0 ** 2)

    def delta_0(self, P, gamma, M_0):
        """
        Equation (2.52b)
        """
        return (P / P_std) * (1 + ((gamma - 1) / 2) * M_0 ** 2) ** (gamma / (gamma - 1))

    def centripetal_force(self, V, R):
        """
        = n

        Equation (2.17)
        """
        n = np.sqrt(1 + (V ** 2 / (g_0 ** 2 * R)) ** 2)
        return n

    def dV_dt(self, V_final, V_initial, delta_t_allowable):
        return ((V_final - V_initial) / delta_t_allowable)

    def s_G(self, alpha, beta, C_Lmax, rho, k_TO, W_TO, S, ksi_TO):
        """
        ground roll distance

        Equation (2.25)

        Look at "note" page 30-31
        """
        pass

    def s_TO(self, alpha_wet, beta, C_Lmax, rho, k_TO, T_SL, W_TO, S):
        """
        takeoff distance

        Equation (2.E1)
        :param alpha_wet:
        :param beta:
        :param C_Lmax:
        :param rho:
        :param k_TO:
        :param T_SL:
        :param W_TO:
        :param S:
        """
        return ((k_TO ** 2 * beta ** 2) /
                (rho * C_Lmax * alpha_wet * (T_SL / W_TO))) * (W_TO / S) + (
                t_R * k_TO) * np.sqrt((2 * beta) / (rho * C_Lmax)) * np.sqrt(W_TO / S)

    def s_R(self, beta, C_Lmax, rho, k_TO, W_TO, S):
        return t_R * k_TO * np.sqrt((2 * beta) / (rho * C_Lmax)) * np.sqrt(W_TO / S)

    def C_D(self, K1, K2, C_L, C_D0):
        """
        Equation (2
        :param K1:
        :param K2:
        :param C_L:     lift coefficient
        :param C_D0:    drag coefficient at zero lift
        :return:
        """
        return K1 * C_L ** 2 + K2 * C_L + C_D0


class MasterEqn(CommonFunctionality):

    def installed_thrust(self, alpha, T_SL):
        """
        = T

        Equation (2.3)
        :param alpha:   installed full throttle thrust lapse
        :param T_SL:    thrust loading at sea level
        """
        return alpha * T_SL

    def instantaneous_weight(self, beta, W_TO):
        """
        = W

        Equation (2.4)
        :param beta:    instantaneous weight fraction
        :param W_TO:    weight at takeoff
        """
        return beta * W_TO

    def master_thrust_to_weight_v1(self, alpha, beta, q, n, K1, K2, C_D0, C_DR, V, dV_dt, dh_dt, W_TO, S):
        """
        = T_SL / W_TO

        Equation (2.11)
        :param q:
        :param n:       load factor, number of g's
        :param K1:
        :param K2:
        :param C_D0:    drag coefficient at zero lift
        :param C_DR:    drag coefficient from external sources
        :param V:       velocity
        :param dV_dt:   acceleration
        :param dh_dt:   rise / descent
        :param W_TO:    weight at takeoff
        :param S:       surface of wing?
        """
        P_s = self.P_s(V, dh_dt, dV_dt)
        return (beta / alpha) * ((q * S) / (beta * W_TO) * (
                K1 * ((n * beta * W_TO) / (q * S)) ** 2
                + K2 * ((n * beta * W_TO) / (q * S))
                + C_D0 + C_DR) + (P_s / V))

    def master_thrust_to_weight_v2(self, alpha, beta, q, n, K1, K2, C_D0, C_DR, V, dV_dt, dh_dt, wing_loading):
        """
        = T_SL / W_TO

        Equation (2.11)
        :param wing_loading:    W_TO / S
        """
        P_s = self.P_s(V, dh_dt, dV_dt)
        return (beta / alpha) * (q / (beta * wing_loading) * (
                K1 * ((n * beta * wing_loading) / q) ** 2
                + K2 * ((n * beta * wing_loading) / q)
                + C_D0 + C_DR) + (P_s / V))


class Case1(MasterEqn):
    """
    Constant Altitude / Speed Cruise (P_s = 0)
    Given: dh/dt=0, dV/dt=0, n=1, h, V, q

    Page 44: Supercruise: C_DR = K2 = 0
    Page 51: Mission phases 6-7 and 8-9: Supersonic penetration and escape dash: C_DR = K2 = 0
    Page 54: Maximum Mach number: C_DR = K2 = 0
    """

    def wing_loading_min(self, beta, q, K1, C_D0, C_DR):
        return (q / beta) * np.sqrt((C_D0 + C_DR) / K1)

    def thrust_to_weight_min(self, alpha, beta, n, K1, K2, C_D0, C_DR):
        """
        Corresponds to maximum range
        """
        return ((n * beta) / alpha) * (
                2 * np.sqrt((C_D0 + C_DR) * K1) + K2)

    def thrust_to_weight(self, alpha, beta, q, K1, K2, C_D0, C_DR, V, wing_loading):
        n = 1
        dh_dt = 0
        dV_dt = 0
        return self.master_thrust_to_weight_v2(alpha, beta, q, n, K1, K2, C_D0, C_DR, V, dV_dt, dh_dt, wing_loading)


class Case2(MasterEqn):
    """
    Constant Speed Climb (P_s = dh/dt)
    Given: dV/dt=0, n=1, h, dh/dt>0, V, q

    same wing_loading_min, thrust_to_weight_min as Case1
    """

    def wing_loading_min(self, beta, q, K1, C_D0, C_DR):
        return (q / beta) * np.sqrt((C_D0 + C_DR) / K1)

    def thrust_to_weight_min(self, alpha, beta, K1, K2, C_D0, C_DR, V, dh_dt):
        dV_dt = 0
        P_s = self.P_s(V, dh_dt, dV_dt)
        return (beta / alpha) * (
                2 * np.sqrt((C_D0 + C_DR) * K1) + K2 + (P_s / V))

    def thrust_to_weight(self, alpha, beta, q, K1, K2, C_D0, C_DR, V, dh_dt, wing_loading):
        n = 1
        dV_dt = 0
        return self.master_thrust_to_weight_v2(alpha, beta, q, n, K1, K2, C_D0, C_DR, V, dV_dt, dh_dt, wing_loading)


class Case3(MasterEqn):
    """
    Constant Altitude / Speed Turn (P_s = 0)
    Given: dh/dt=0, dV/dt=0, h, V, q, n>1

    same wing_loading_min, thrust_to_weight_min as Case1
    """

    def n(self, R):
        """"
        :param R: radius of turn
        """
        return self.centripetal_force()

    def thrust_to_weight(self, alpha, beta, q, n, K1, K2, C_D0, C_DR, V, wing_loading):
        dV_dt = 0
        dh_dt = 0
        return self.master_thrust_to_weight_v2(alpha, beta, q, n, K1, K2, C_D0, C_DR, V, dV_dt, dh_dt, wing_loading)


class Case4(MasterEqn):
    """
    Horizontal Acceleration [P_s = (V/g0)(dV/dt)]
    Given: dh/dt=0, n=1, h, V_initial, V_final, delta_t_allowable

    Page 52: Mission phase 7-8: Horizontal acceleration: C_DR = K2 = 0
    """

    # def delta_t_function(self, V):
    #     return (1 / (2 * g_0)) * (V / self.P_s)
    #
    # def delta_t(self,):
    #     result, error = integrate.quad(delta_t_function, self.V_initial, self.V_final)
    #     delta_t_allowable = result
    #     return delta_t_allowable

    def dV_dt(self, V_final, V_initial, delta_t_allowable):
        """
        :param V_final: final velocity
        :param V_initial: initial velocity
        :param delta_t_allowable: allowable time to accelerate
        """
        return self.dV_dt(V_final, V_initial, delta_t_allowable)

    def thrust_to_weight(self, alpha, beta, q, K1, K2, C_D0, C_DR, V, dV_dt, wing_loading):
        n = 1
        dh_dt = 0
        return self.master_thrust_to_weight_v2(alpha, beta, q, n, K1, K2, C_D0, C_DR, V, dV_dt, dh_dt, wing_loading)

    # ????????????????
    # ????????????????
    # ????????????????


class Case5(MasterEqn):
    """
    Takeoff Ground Roll (s_G), when T_SL >> (D + R)
    Given: dh/dt=0, s_G, rho, C_Lmax, V_TO = k_TO V_STALL
    """

    def thrust_to_weight(self, alpha, beta, C_Lmax, rho, s_G, k_TO, wing_loading):
        """
        Equation (2.22)
        :param C_Lmax:  maximum coefficient of lift
        :param rho:     density?
        :param s_G:     ground roll distance
        :param k_TO:
        """
        return ((beta ** 2 / alpha) *
                (k_TO ** 2 / (s_G * rho * g_0 * C_Lmax)) *
                wing_loading)

    def wing_loading(self, alpha_wet, beta, C_Lmax, rho, k_TO, T_SL, W_TO):
        """
        Equation (2.E2)
        """
        # a = (k_TO ** 2 * beta ** 2) / (rho * C_Lmax * alpha_wet * (thrust_to_weight(wing_loading)))
        # b = (t_R * k_TO) * np.sqrt((2 * beta) / (rho * C_Lmax))
        # s_TO = ((k_TO ** 2 * beta ** 2) /
        #         (rho * C_Lmax * alpha_wet * (thrust_to_weight(wing_loading)))) * (wing_loading) + (
        #         t_R * k_TO) * np.sqrt((2 * beta) / (rho * C_Lmax)) * np.sqrt(wing_loading)
        # c = s_TO
        a = (k_TO ** 2 * beta ** 2) / (rho * C_Lmax * alpha_wet * (T_SL / W_TO))
        b = (t_R * k_TO) * np.sqrt((2 * beta) / (rho * C_Lmax))
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

    """
    :param C_D:     drag coefficient
    :param C_L:     lift coefficient
    :param mu_TO:
    """

    # self.R = q * C_DR * S + mu_TO * (beta * W_TO - q * C_L * S)
    # self.ksi_TO = C_D + C_DR - mu_TO * C_L

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

    # self.D = q * C_D
    # self.R = q * C_DR * S + mu_TO * (beta * W_TO - q * C_L * S)
    # self.V_TD = k_TD * V_STALL

    def case7_thrust_to_weight(self, alpha, beta, C_Lmax, rho, s_B, k_TD, W_TO, S):
        """
        Equation (2.34)
        Reverse thrust
        """
        return ((beta ** 2 / (- alpha)) *
                (k_TD ** 2 / (s_B * rho * g_0 * C_Lmax)) *
                (W_TO / S))
    # ????????????????
    # ????????????????
    # ????????????????


class Case8(MasterEqn):
    """
    Service Ceiling (P_s = dh/dt)
    Given: dV/dt=0, n=1, h, dh/dt>0, C_L
    """

    def thrust_to_weight(self, alpha, beta, K1, K2, C_D0, C_DR, V, dh_dt, C_L):
        return ((beta / alpha) *
                (K1 * C_L + K2 + ((C_D0 + C_DR) / C_L) +
                 (1 / V) * dh_dt))


class Case9(MasterEqn):
    """
    Takeoff Climb Angle
    Given: theta, n=1, dV/dt=0, C_DR, C_Lmax, k_TO, h, sigma
    """

    def thrust_to_weight(self, alpha, beta, K1, K2, C_D0, C_DR, C_Lmax, k_TO, theta):
        """
        Equation (2.43)
        :param theta:
        """
        return ((beta / alpha) *
                ((K1 * C_Lmax / k_TO ** 2) + K2 +
                 ((C_D0 + C_DR) / (C_Lmax / k_TO ** 2)) + np.sin(theta)))

    def V(self, beta, C_Lmax, k_TO, sigma, W_TO, S):
        """
        Equation (2.44)
        :param sigma:
        """
        return np.sqrt(((2 * beta * k_TO ** 2) / (sigma * rho_SL * C_Lmax)) * (W_TO / S))


class Mach_vs_ThrustLapse(CommonFunctionality):
    """
    Functionality to create plot on page 42
    """

    def __init__(self, h, gamma, M_0, TR):
        """
        :param h:       altitude
        :param M_0:     mach number
        :param TR:      throttle ratio
        """
        self.h = h
        self.gamma = gamma
        self.T = self.temperature(h)
        self.P = self.pressure(h)
        self.M_0 = M_0
        self.TR = TR
        self.theta_0 = self.theta_0(self.T, self.gamma, self.M_0)
        self.delta_0 = self.delta_0(self.P, self.gamma, self.M_0)

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
