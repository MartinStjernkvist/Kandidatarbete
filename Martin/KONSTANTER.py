"""in m"""
r_tip_fan = 0.9
r_hub_fan = r_tip_fan - 0.62
r_splitter = r_tip_fan - 0.316
r_hub_splitter = r_hub_fan

M_cruise = 1.7
h_cruise = 18288
gamma_a = 1.4
gamma_c = 1.33
g_0 = 9.81  # gravitation
rho_SL = 1.225  # kg / m^3 air sealevel
t_R = 3  # runtime on takeoff
s_TO = 1500  # takeoff distance in ft

T_std = 288.15  # K
P_std = 101_325  # Pa
R = 287.05  # J/kg*K

T_turbine = 1650  # K
gamma_t = ...
g_c = 1
A_fan = 2.296  # m^2

F_symphony = 160_000  # N

P_effektivitet = 0.99  # Pt2 Verklig / Pt2 teoretisk

# s_G = ? # ground roll distance
# s_B = ? # breaking distance
# s_TO = ? # takeoff distance