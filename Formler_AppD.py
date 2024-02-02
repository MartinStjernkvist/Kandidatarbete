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
placeholder = 10 ** 20  # use as placeholder when creating instances of classes,


# not used in calculations of the instance

class AppendixD():

    def __init__(self, T_t_4_max, T_t_SLS, gamma_c, eta_c, eta_m, beta, c_pt, c_pc, T_0, M_0, T_t_4, f, tau_c,
                 M_0_break):
        self.T_t_4_max = T_t_4_max
        self.T_t_SLS = T_t_SLS
        self.gamma_c = gamma_c
        self.eta_c = eta_c
        self.eta_m = eta_m
        self.c_pt = c_pt
        self.c_pc = c_pc
        self.M_0 = M_0
        self.T_0 = T_0
        self.beta = beta
        self.f = f
        self.T_t_4 = T_t_4
        self.tau_c = tau_c
        self.M_0_break = M_0_break

    def theta_0_break(self, T_t_4_max, T_t_SLS):
        return T_t_4_max / T_t_SLS

    def theta_0(self, T_0, gamma_c, M_0):
        """
        Equation (D.1)
        Same as (2.52a), but with gamma_c
        """
        theta_tau_r = (T_0 / T_std) * (1 + ((gamma_c - 1) / 2) * M_0 ** 2)
        return theta_tau_r

    def tau_c(self, eta_m, beta, f, c_pt, T_t_4, c_pc, theta_0):
        """
        Compressor total temperature ratio
        Equation (D.2)

        Depends only on:
        - turbine total temperature ratio tau_t
        - throttle setting: T_t_4
        - flight conditions theta_0
        """
        tau_c = (1 + eta_m * (1 - beta) *
                 ((1 + f) * (1 / T_std) *
                  ((c_pt * T_t_4) / (c_pc * theta_0))))
        return tau_c

    def ThrottleRatio(self, T_t_4_max, T_t_SLS):
        """
        Equation (D.6)
        """
        theta_0_break = T_t_4_max / T_t_SLS
        return theta_0_break

    def MachBreak(self, gamma_c, theta_0_break):
        """
        Equation (D.7)
        :param gamma_c:
        :param theta_0_break:
        :return:
        """
        M_0_break = np.sqrt((2 / (gamma_c - 1)) * (theta_0_break - 1))
        return M_0_break


class AppendixD_dup(AppendixD):
    def __init__(self, T_t_4_max, T_t_SLS, gamma_c, eta_c, eta_m, beta, c_pt, c_pc, T_0, M_0, T_t_4, f, tau_c,
                 M_0_break):
        super().__init__(T_t_4_max, T_t_SLS, gamma_c, eta_c, eta_m, beta, c_pt, c_pc, T_0, M_0, T_t_4, f, tau_c,
                         M_0_break)
        self.theta_0_break = self.theta_0_break(T_t_4_max, T_t_SLS)
        self.M_0_break = self.MachBreak(gamma_c, self.theta_0_break)
        self.theta_0 = self.theta_0(T_0, gamma_c, M_0)
        self.tau_c = self.tau_c(eta_m, beta, f, c_pt, T_t_4, c_pc, self.theta_0)




    # def calc_theta_0_break_dup(self, T_t_4_max, T_t_SLS):
    # def calc_theta_0_break_dup(self, T_t_4_max, T_t_SLS):
    #     return T_t_4_max / T_t_SLS


    def get_T_t_4_max(self):
        return self.T_t_4_max

    def get_T_t_SLS(self):
        return self.T_t_SLS

    def get_gamma_c(self):
        return self.gamma_c

    def get_eta_c(self):
        return self.eta_c

    def get_eta_m(self):
        return self.eta_m

    def get_beta(self):
        return self.beta

    def get_c_pt(self):
        return self.c_pt

    def get_c_pc(self):
        return self.c_pc

    def get_T_0(self):
        return self.T_0

    def get_M_0(self):
        return self.M_0

    def get_T_t_4(self):
        return self.T_t_4

    def get_f(self):
        return self.f

    def get_tau_c(self):
        return self.tau_c

    def get_M_0_break(self):
        return self.M_0_break

    def get_theta_0(self):
        return self.theta_0_break

    # def theta_0_dup(self):
    #     """
    #     Equation (D.1)
    #     Same as (2.52a), but with gamma_c
    #     """
    #     theta_tau_r = (self.T_0 / T_std) * (1 + ((self.gamma_c - 1) / 2) * self.M_0 ** 2)
    #     return theta_tau_r
    #
    # def CompressorTotalTemperatureRatioToYield_dup(self):
    #     """
    #     Equation (D.2)
    #
    #     Depends only on:
    #     - turbine total temperature ratio tau_t
    #     - throttle setting: T_t_4
    #     - flight conditions theta_0
    #     """
    #     tau_c = 1 + self.eta_m * (1 - self.beta) * ((1 + self.f) * (1 / T_std) *
    #                                                 ((self.c_pt * self.T_t_4) / (self.c_pc * self.theta_0)))
    #
    #     return tau_c
    #
    # def ThrottleRatio_dup(self):
    #     """
    #     Equation (D.6)
    #     """
    #     theta_0_break = self.T_t4max / self.T_tSLS
    #     return theta_0_break
    #
    # def MachBreak_dup(self):
    #     """
    #     Equation (D.7)
    #     """
    #     M_0_break = np.sqrt((2 / (self.gamma_c - 1)) * (self.theta_0_break - 1))
    #     return M_0_break


T_t_4_max = 1
T_t_SLS = 1
gamma_c = 1.5
eta_c = 1
eta_m = 1
beta = 1
c_pt = 1
c_pc = 1
T_0 = 1
M_0 = 1
T_t_4 = 1
f = 1
tau_c = 1
M_0_break = 1

inst = AppendixD_dup(T_t_4_max, T_t_SLS, gamma_c, eta_c, eta_m, beta, c_pt, c_pc, T_0, M_0, T_t_4, f, tau_c, M_0_break)

new_T_t_4_max = inst.get_T_t_4_max()
new_T_t_SLS = inst.get_T_t_SLS()
new_gamma_c = inst.get_gamma_c()
new_eta_c = inst.get_eta_c()
new_eta_m = inst.get_eta_m()
new_beta = inst.get_beta()
new_c_pt = inst.get_c_pt()
new_c_pc = inst.get_c_pc()
new_T_0 = inst.get_T_0()
new_M_0 = inst.get_M_0()
new_T_t_4 = inst.get_T_t_4()
new_f = inst.get_f()
new_tau_c = inst.get_tau_c()
new_M_0_break = inst.get_M_0_break()

print(
    f'\nT_t_4_max {T_t_4_max}, \nT_t_SLS {T_t_SLS}, \ngamma_c {gamma_c}, \neta_c {eta_c}, \neta_m {eta_m}, \nbeta {beta}, \nc_pt {c_pt}, \nc_pc {c_pc}, \nT_0 {T_0}, \nM_0 {M_0}, \nT_t_4 {T_t_4}, \nf {f}, \ntau_c {tau_c}, \nM_0_break {M_0_break}')
print(
    f'\nnew_T_t_4_max {new_T_t_4_max}, \nnew_T_t_SLS {new_T_t_SLS}, \nnew_gamma_c {new_gamma_c}, \nnew_eta_c {new_eta_c}, \nnew_eta_m {new_eta_m}, \nnew_beta {new_beta}, \nnew_c_pt {new_c_pt}, \nnew_c_pc {new_c_pc}, \nnew_T_0 {new_T_0}, \nnew_M_0 {new_M_0}, \nnew_T_t_4 {new_T_t_4}, \nnew_f {new_f}, \nnew_tau_c {new_tau_c}, \nnew_M_0_break {new_M_0_break}')
