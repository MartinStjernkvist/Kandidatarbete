import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
from FORM_Exempelupg import Ex
from FORM_8Performance import Perf
from FORM_General import Gen
from FORM_DesignTask2 import DT2
from FORM_Exempelupg import Ex
from KONSTANTER import *
from scipy.optimize import fsolve
"""
Lite syntax

P = tryck
T = temperatur
a = ljudhastighet
C_p = värmekapacitet (konstant tryck)
C_v = värmekapacitet (konstant temperatur)
gamma = specifik värme ratio
omega = vinkelhastighet
U = bladhastighet
c = lufthastighet
MFP = mass flow parameter
psi = steglast
Delta_H = entalpiändring
Delta_T = temperaturändring
"""

class Stage():
    def __init__(self, T_t, M_ax, M_rel, r_tip_list, r_hub_list, stage_loading, efficiency, ):
        self.T_t = T_t
        self.M_ax = M_ax
        self.M_rel = M_rel
        self.r_tip_list = r_tip_list
        self.r_hub_list = r_hub_list
        self.stage_loading = stage_loading
        self.efficiency = efficiency


        self.C_p = Perf().C_p_air(self.T_t)
        self.gamma = Perf().gamma(self.C_p, R)
        self.A = Gen().A(self.r_tip_list[1], self.r_hub_list[1])
        self.MFP = Gen().MFP(self.M)



