import numpy as np

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

class AppendixD:

    def theta_0(self, T_0, gamma_c, M_0):
        """
        dimensionless ratio of freestream total temperature to
        sea level static temperature of the standard atmosphere

        = theta_tau_r

        Equation (D.1)
        Same as (2.52a), but with gamma_c

        :param T_0:         freestream temperature
        :param gamma_c:
        :param M_0:
        """
        theta_tau_r = (T_0 / T_std) * (1 + ((gamma_c - 1) / 2) * M_0 ** 2)
        return theta_tau_r

    def tau_c(self, eta_m, beta, f, c_pt, T_t_4, c_pc, theta_0):
        """
        compressor total temperature ratio

        Equation (D.2)

        Depends only on:
        - turbine total temperature ratio tau_t
        - throttle setting: T_t_4
        - flight conditions theta_0

        :param eta_m:
        :param beta:
        :param f:
        :param c_pt:
        :param T_t_4:
        :param c_pc:
        :param theta_0:
        """
        return (1 + eta_m * (1 - beta) *
                ((1 + f) * (1 / T_std) *
                 ((c_pt * T_t_4) / (c_pc * theta_0))))

    def theta_0_break(self, T_t_4_max, T_t_SLS):
        """
        throttle ratio

        Equation (D.6)

        :param T_t_4_max:
        :param T_t_SLS:
        """
        return T_t_4_max / T_t_SLS

    def M_0_break(self, gamma_c, theta_0_break):
        """
        mach break

        Equation (D.7)

        :param gamma_c:
        :param theta_0_break:
        """
        return np.sqrt((2 / (gamma_c - 1)) * (theta_0_break - 1))

"""
T_t_4_max = 
T_t_SLS = 
gamma_c = 
eta_c = 
eta_m = 
beta = 
c_pt = 
c_pc = 
T_0 = 
M_0 = 
T_t_4 = 
f = 
tau_c = 
M_0_break = 
"""

theta_0_break = AppendixD().theta_0_break(1, 1)
print(theta_0_break)

# print(
#     f'\nT_t_4_max {T_t_4_max}, \nT_t_SLS {T_t_SLS}, \ngamma_c {gamma_c}, \neta_c {eta_c}, \neta_m {eta_m}, \nbeta {beta}, \nc_pt {c_pt}, \nc_pc {c_pc}, \nT_0 {T_0}, \nM_0 {M_0}, \nT_t_4 {T_t_4}, \nf {f}, \ntau_c {tau_c}, \nM_0_break {M_0_break}')
# print(
#     f'\nnew_T_t_4_max {new_T_t_4_max}, \nnew_T_t_SLS {new_T_t_SLS}, \nnew_gamma_c {new_gamma_c}, \nnew_eta_c {new_eta_c}, \nnew_eta_m {new_eta_m}, \nnew_beta {new_beta}, \nnew_c_pt {new_c_pt}, \nnew_c_pc {new_c_pc}, \nnew_T_0 {new_T_0}, \nnew_M_0 {new_M_0}, \nnew_T_t_4 {new_T_t_4}, \nnew_f {new_f}, \nnew_tau_c {new_tau_c}, \nnew_M_0_break {new_M_0_break}')
