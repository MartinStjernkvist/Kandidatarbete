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

    def __init__(self, T_t_4_max, T_t_SLS, gamma_c, eta_c, eta_m, beta, c_pt, c_pc, T_0, M_0, T_t_4, f):
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

    def calc_theta_0_break(self, T_t_4_max, T_t_SLS):
        return T_t_4_max / T_t_SLS


    def theta_0(self, T_0, gamma_c, M_0):
        """
        Equation (D.1)
        Same as (2.52a), but with gamma_c
        """
        theta_tau_r = (T_0 / T_std) * (1 + ((gamma_c - 1) / 2) * M_0 ** 2)
        return theta_tau_r

    def CompressorTotalTemperatureRatioToYield(self, eta_m, beta, f, c_pt, T_t_4, c_pc, theta_0):
        """
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

    def ThrottleRatio(self):
        """
        Equation (D.6)
        """
        theta_0_break = self.T_t_4_max / self.T_t_SLS
        return theta_0_break

    def MachBreak(self):
        """
        Equation (D.7)
        """
        M_0_break = np.sqrt((2 / (self.gamma_c - 1)) * (self.theta_0_break - 1))
        return M_0_break


class AppendixD_dup(AppendixD):
    def __init__(self, T_t_4_max, T_t_SLS, gamma_c, eta_c, eta_m, beta, c_pt, c_pc, T_0, M_0, T_t_4, f):
        super().__init__(T_t_4_max, T_t_SLS, gamma_c, eta_c, eta_m, beta, c_pt, c_pc, T_0, M_0, T_t_4, f)
        self.theta_0_break = self.calc_theta_0_break(T_t_4_max, T_t_SLS)
        self.theta_0 = self.theta_0(T_0, gamma_c, M_0)
        self.tau_c = self.CompressorTotalTemperatureRatioToYield(eta_m, beta, f, c_pt, T_t_4, c_pc, self.theta_0)

    # def calc_theta_0_break_dup(self, T_t_4_max, T_t_SLS):
    # def calc_theta_0_break_dup(self, T_t_4_max, T_t_SLS):
    #     return T_t_4_max / T_t_SLS

    def get_theta_0_dup(self):
        return self.theta_0_break

    def get_tau_c(self):
        return self.tau_c

    def theta_0_dup(self):
        """
        Equation (D.1)
        Same as (2.52a), but with gamma_c
        """
        theta_tau_r = (self.T_0 / T_std) * (1 + ((self.gamma_c - 1) / 2) * self.M_0 ** 2)
        return theta_tau_r

    def CompressorTotalTemperatureRatioToYield_dup(self):
        """
        Equation (D.2)

        Depends only on:
        - turbine total temperature ratio tau_t
        - throttle setting: T_t_4
        - flight conditions theta_0
        """
        tau_c = 1 + self.eta_m * (1 - self.beta) * ((1 + self.f) * (1 / T_std) *
                                                    ((self.c_pt * self.T_t_4) / (self.c_pc * self.theta_0)))

        return tau_c

    def ThrottleRatio_dup(self):
        """
        Equation (D.6)
        """
        theta_0_break = self.T_t4max / self.T_tSLS
        return theta_0_break

    def MachBreak_dup(self):
        """
        Equation (D.7)
        """
        M_0_break = np.sqrt((2 / (self.gamma_c - 1)) * (self.theta_0_break - 1))
        return M_0_break


inst = AppendixD_dup(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
theta_0 = inst.get_theta_0_dup()
tau_c = inst.get_tau_c()
print(theta_0, tau_c)

new_tau_c = inst.CompressorTotalTemperatureRatioToYield_dup()
print(new_tau_c)
