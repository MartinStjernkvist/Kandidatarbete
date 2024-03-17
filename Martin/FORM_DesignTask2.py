import numpy as np
from KONSTANTER import *


class DT2():
    def M_ax(self, M, theta):
        """
        axial mach number, same as M_cruise at cruise
        :param theta:   angle of attack
        """
        return M * np.cos(theta)

    def hub_tip_ratio(self, r_hub, r_tip):
        """
        hub tip ratio
        """
        return r_hub / r_tip

    def AR(self, r_hub_1, r_hub_2, r_tip_1, r_tip_2):
        """
        aspect ratio
        """
        return ((r_tip_1 + r_tip_2) / 2 - (r_hub_1 + r_hub_2) / 2) / 2

    def l_ax(self, alpha, r_tip_1, r_hub_1, AR):
        return (1.98 * r_tip_1 - 2 * r_hub_1) / (2 * AR + np.tan(alpha))

    def U_tip(self):
        pass

    def A_duct_entry(self, A_2, BPR):
        return A_2 / (BPR + 1)

    def U(self, M_rel, a, c_air):
        """
        blade speed
        """
        return np.sqrt((M_rel * a)**2 - c_air**2)

    def U_mid(self, r_tip, r_hub, omega):
        return ((r_tip + r_hub) / 2) * omega

    # def stage_loading(self, Delta_H, U_mid_sum):
    #     return 2 * Delta_H / U_mid_sum

    def delta_H(self, psi, sum_U_mid_squared):
        return psi * sum_U_mid_squared / 2
