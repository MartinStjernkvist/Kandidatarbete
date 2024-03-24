import numpy as np
from scipy.optimize import fsolve
from KONSTANTER import *


class Gen:
    """
    General equations, mostly from Mattingly
    """

    def P_s(self, V, dh_dt, dV_dt):
        """
        Weight specific excess power

        Mattingly (2.2b)
        """
        return dh_dt + (V / g_0) * dV_dt

    def temperature(self, h):
        """
        Introkompendium (1.29a)
        """
        if h <= 11000:
            T_0 = T_std - (6.5 * 10 ** (-3) * h)
        else:
            T_0 = T_std - (6.5 * 10 ** (-3) * 11000)
        return T_0

    def pressure(self, h):
        """
        Introkompendium (1.29b)
        """
        # if h <= 11000:
        #     P_0 = P_std * ((T_std - 6.5 * 10 ** (-3) * h) / T_std) ** (g_0 / (6.5 * 10 ** (-3) * R))
        #     return P_0
        # else:
        #     P_0 = P_std * ((T_std - 6.5 * 10 ** (-3) * 11000) / T_std) ** (g_0 / (6.5 * 10 ** (-3) * R))
        #     return P_0
        if h <= 11000:
            P_0 = P_std * ((T_std - 6.5 * 10 ** (-3) * h) / T_std) ** (g_0 / (6.5 * 10 ** (-3) * R))
            return P_0
        else:
            return f'altitude exceeds 11km'

    def theta_0(self, T, gamma, M_0):
        """
        dimensionless ratio of freestream total temperature to
        sea level static temperature of the standard atmosphere

        Mattingly (2.52a)
        """
        return (T / T_std) * (1 + ((gamma - 1) / 2) * M_0 ** 2)

    def delta_0(self, P, gamma, M_0):
        """
        Mattingly (2.52b)
        """
        return (P / P_std) * (1 + ((gamma - 1) / 2) * M_0 ** 2) ** (gamma / (gamma - 1))

    def centripetal_force(self, V, R):
        """
        = n

        Mattingly (2.17)
        """
        n = np.sqrt(1 + (V ** 2 / (g_0 ** 2 * R)) ** 2)
        return n

    def dV_dt(self, V_final, V_initial, delta_t_allowable):
        return ((V_final - V_initial) / delta_t_allowable)

    def s_G(self, alpha, beta, C_Lmax, rho, k_TO, W_TO, S, ksi_TO):
        """
        ground roll distance

        Mattingly (2.25)

        Look at "note" page 30-31
        """
        pass

    def s_TO(self, alpha_wet, beta, C_Lmax, rho, k_TO, T_SL, W_TO, S):
        """
        takeoff distance

        Mattingly (2.E1)
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
        Mattingly (2
        :param K1:
        :param K2:
        :param C_L:     lift coefficient
        :param C_D0:    drag coefficient at zero lift
        :return:
        """
        return K1 * C_L ** 2 + K2 * C_L + C_D0

    # def a(self, T, gamma):
    #     """
    #     speed of sound
    #     """
    #     speed_of_sound = np.sqrt(gamma * R * T)
    #     return speed_of_sound

    # def stagnation_temperature(self, T, M, gamma):
    #     return T * (1 + (gamma - 1) / 2 * M ** 2)
    #
    # def stagnation_pressure(self, P, M, gamma):
    #     return P * (1 + (gamma - 1) / 2 * M ** 2) ** (gamma / (gamma - 1))

    def MFP(self, M, gamma, g_c=1, ):
        """
        mass flow parameter

        Mattingly (1.3)
        :param M:       mach number
        :param g_c:     gas constant?
        :param gamma:
        """
        return M * np.sqrt((gamma * g_c / R)) * (1 + ((gamma - 1) / 2) * M ** 2) ** ((gamma + 1) / (2 * (1 - gamma)))

    def mass_flow(self, MFP, P, T, A):
        """
        mdot

        Mattingly (1.3)
        """
        return MFP * P * A / np.sqrt(T)

    def M_solve(self, gamma, MFP):
        def equation_to_solve(M):
            return self.MFP(M, gamma, R) - MFP

        initial_guess = 0
        solution = fsolve(equation_to_solve, initial_guess)
        return solution

    def RPM(self, omega):
        return omega * 60 / (2 * np.pi)

    def mass_flow_bypass(self, mass_flow, BPR):
        return mass_flow * BPR / (1 + BPR)

    def mass_flow_core(self, mass_flow, mass_flow_bypass):
        return mass_flow - mass_flow_bypass

    def A(self, r_tip, r_hub):
        return np.pi * r_tip ** 2 * (1 - r_hub / r_tip)
