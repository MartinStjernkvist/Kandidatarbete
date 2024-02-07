from Kandidatarbete.Martin.FORMLER_Exempel import Ex
# from Thrust import F

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

h = 10668
M_0 = 0.78

T = Ex().temperature(h)
p = Ex().pressure(h)
a = Ex().a(gamma, R, T)
c = M_0 * a
rho = Ex().rho(p, R, T)
print(c)

mdot = 185

A = mdot / (rho * c)
print(A)


overall_PR = 45
fan_PR = 1.48
bypass_PR = 12.5
HPC_PR = 10
turbine_INT = 1650  # K
fan_IPC_HPC_poly = 0.9
HPT_poly = 0.85
LPT_poly = 0.9
combustor_PL = 0.04
tot_AMF = 185  # kg/s
A_0 = 2.1081269749763853  # from M_0 = 0.78, mass flow = 185 kg/s, h = 10668

IPC_PR = overall_PR / (fan_PR * HPC_PR)

print(fan_PR * HPC_PR * IPC_PR)