import numpy as np
from KONSTANTER import *


class DT2():
    def M_ax(self, M, theta):
        """
        axial mach number, same as M_cruise at cruise
        :param theta:   angle of attack
        """
        return M * np.cos(theta)

    def nu(self, r_hub, r_tip):
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

    def U_mid(self, r_tip, r_hub):
        return (r_tip + r_hub) / 2

    def stage_loading(self, Delta_H, U_mid_list):
        return 2 * Delta_H / np.sum([U_mid_list[i] ** 2 for i in range(U_mid_list)])
