import numpy as np
from FORM_Mattingly import Gen
from KONSTANTER import *


class Ex:

    def thrust(self, mdot_0, mdot_f, p_0, p_8, c_0, c_8, A_8):
        """
        F

        Equation (1.5)
        """
        return (mdot_0 + mdot_f) * c_8 + (p_8 - p_0) * A_8 - mdot_0 * c_0

    def rho(self, p, R, T):
        """
        density

        Equation (1.8)
        """
        return p / (R * T)

    def mdot(self, rho, A, c):
        """
        mass flow

        Equation (1.18)
        """
        return rho * A * c

    def a(self, gamma, R, T):
        """
        speed of sound

        Equation (1.22)
        """
        return np.sqrt(gamma * R * T)

    def M(self, a, c):
        """
        Equation (1.23)
        """
        return c / a

    def stagnation_PR(self, gamma, M):
        """
        stagnation pressure ratio
        Equation (1.28b)
        """
        return (1 + ((gamma - 1) / 2) * M ** 2) ** (gamma / (gamma - 1))

    def stagnation_PR_fcn_TR(self, T_final, T_initial, eta, gamma):
        temperature_ratio = T_final / T_initial
        return temperature_ratio ** (gamma * eta / (gamma - 1))

    def stagnation_TR(self, gamma, M):
        """
        stagnation temperature ratio
        """
        T_0_T = 1 + ((gamma - 1) / 2) * M ** 2
        return T_0_T

    def stagnation_TR_compressor(self, gamma, M, eta_pc):
        """
        stagnation temperature ratio compressor

        Equation (1.32)
        """
        T_03_T_02 = self.stagnation_PR(gamma, M) ** ((gamma - 1) / (gamma * eta_pc))
        return T_03_T_02

    def stagnation_TR_turbine(self, gamma, M, eta_pt):
        """
        stagnation temperature ratio turbine
        Equation (1.35)
        """
        T_04_T_05 = self.stagnation_PR(gamma, M) ** ((gamma - 1) * eta_pt / gamma)
        return T_04_T_05

    def p_0(self, p, gamma, M):
        return p * (1 + ((gamma - 1) / 2) * M ** 2) ** (gamma / (gamma - 1))

    def p(self, p_0, gamma, M):
        return p_0 / (1 + ((gamma - 1) / 2) * M ** 2) ** (gamma / (gamma - 1))

    def T_0(self, T, gamma, M):
        return T * (1 + ((gamma - 1) / 2) * M ** 2)

    def T(self, T_0, gamma, M):
        return T_0 / (1 + ((gamma - 1) / 2) * M ** 2)

    def W_compressor(self, mdot, c_p, T_03, T_02):
        """
        work compressor
        negative value

        Equation (1.33)
        """
        return - mdot * c_p * (T_03 - T_02)

    def W_turbine(self, mdot, mdot_f, c_p, T_04, T_05):
        """
        work turbine

        Equation (1.37)
        """
        return (mdot + mdot_f) * c_p * (T_04 - T_05)

    def SFC(self, mdot_f, F):
        """
        specific fuel consumption
        """
        return mdot_f / F

    def useful_power(self, F, c_0):
        """
        useful power
        :param c_0:     speed
        """
        P = F * c_0
        return P

    def eta_0(self, F, c_0, mdot_f, LHV):
        return F * c_0 / (mdot_f * LHV)

    def eta_thermal(self, Wdot_kin, mdot_f, LHV):
        """
        thermal efficiency
        :param LHV:     lower heating value
        """
        return Wdot_kin / (mdot_f * LHV)

    def eta_propulsive(self, Wdot_kin, F, c_0):
        """
        propulsive efficiency
        """
        return F * c_0 / Wdot_kin

    def temperature(self, h):
        """
        Introkompendium: Equation (1.29a)
        """
        if h <= 11000:
            T_0 = T_std - (6.5 * 10 ** (-3) * h)
        else:
            T_0 = T_std - (6.5 * 10 ** (-3) * 11000)
        return T_0

    def pressure(self, h):
        """
        Introkompendium: Equation (1.29b)
        """
        # if h <= 11000:
        #     P_0 = P_std * ((T_std - 6.5 * 10 ** (-3) * h) / T_std) ** (g_0 / (6.5 * 10 ** (-3) * R))
        #     return P_0
        # else:
        #     P_0 = P_std * ((T_std - 6.5 * 10 ** (-3) * 11000) / T_std) ** (g_0 / (6.5 * 10 ** (-3) * R))
        #     return P_0
        P_0 = P_std * ((T_std - 6.5 * 10 ** (-3) * h) / T_std) ** (g_0 / (6.5 * 10 ** (-3) * R))
        return P_0
