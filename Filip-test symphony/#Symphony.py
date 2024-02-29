#Symphony
#-----------------------------------Almänna värden-------------------------------------------------------
import numpy as np

T_SLS = 288.15 #kelvin Sea level standar
P_SLS = 101325 #Pa Sea level standard
g = 9.81 #Jordens gravitation
R = 287.05 # luftens unika R
gamma_a = 1.4 # antag att gamam för luft är 1.4
cp_a = gamma_a * R / (gamma_a - 1)  # Obs! Detta gäller för ideala gaser
M_0 = 1.7 #Detta är cruise punkt
M_2 = 0.6 #Detta var givet enligt en källa
n_inlet = 0.98 #Detta beskriver tryckförlust från P_t1 till P_t2
A_2 = 2.296 # Beräknat i parameter värden dokumentet

psi_fläkt = 0.6 #** %från handledning %%%%%%%%%%HITTA VÄRDE%%%%%%%%%%%% Detta är stage-load
M_rel_fläkt = 1.5 #** %från handledning %%%%%%%%%%%%%%HITTA VÄRDE%%%%%%%%%%%%%
n_fläkt = 0.98 # effektivitet för fläkten
r_tipp_fläkt = 0.915 #approximerat, [m]
r_hub_fläkt = 0.29625 #approximerat, [m]


#%%-------------------------------------Beräkningar---------------------------
T_0, P_0, a_0 = Ambient_air(18288) #Beräknar luftegenskaper för omgivningen

T_t2, Pt_2, T_2, P_2 ,mdot, a_2 = Station_2(T_0, P_0)

T_t2i, P_t2i, FPR = Station_2i(T_t2, P_t2, a_2)


#---------------------------------Funktioner--------------------------------------
def Ambient_air(h, T_std = T_SLS, P_std = P_SLS, gamma = gamma_a, R = R):
    """
    Ger tillbaka T_0 och P_0 från T och P vid havsnivå
    """
    h = min(h,11000) #formlerna 1.29 från introkompendie gäller endast till 11000 meter
    T_0 = T_std -6.5*h/(10**3)
    P_0 = P_std * ((T_std - (6.5 * 10 ** (-3) * h)) / T_std) ** (g / (6.5 * 10 ** (-3) * R))
    a_0 = np.sqrt(gamma_a*R*T_0)
    return [T_0, P_0, a_0]

def Station_2(T_0, P_0, A_2 = A_2, M_0 = M_0, M_2 = M_2, gamma = gamma_a, n = n_inlet, R = R): #Innan fläkten
    """
    Ger tillbaka termodynamiska egenskaper innan fläkten + det toala luft-massflödet
    """
    temp = 1 + (gamma - 1) / 2 * M_0 ** 2
    temp2 = 1 + (gamma - 1) / 2 * M_2 ** 2
    T_t2 = T_0 * temp
    P_t2 = P_0 * temp ** (gamma / (gamma-1)) * n
    MFP = M_2 * np.sqrt(gamma / R) * temp2 ** ((gamma + 1) / (2 * (1 - gamma)))  # 1.3 mattingly
    mdot = A_2 * (P_t2 * MFP) / np.sqrt(T_t2)
    T_2 = T_t2 / temp2
    a_2 = np.sqrt(gamma * R * T_2)
    return [T_t2, Pt_2, T_2, P_2, mdot, a_2]

def Station_2i(T_t2, Pt_2, a_2, cp = cp_a, M_2 = M_2, M_rel = M_rel_fläkt, r_tipp = r_tipp_fläkt, r_hub = r_hub_fläkt, psi = psi_fläkt, n = n_fläkt): # Efter fläkten 
    c_2 = M_2 * a_2 
    u_tipp = np.sqrt((M_rel*a_2) ** 2 - (c_2 ** 2)) #bladhastighet vid tippen
    omega = u_tipp / r_tipp; #vinkelhastighet för fläkt
    rpm = ((omega * 60) / (2 * pi)) #rpm för fläkt
    r_mitt = (r_tipp + r_hub) / 2 #radie till berkäning av hastighet u_mitt
    u_mitt = r_mitt * omega #medel av hastighet utöver fläkten 
    delta_h = (psi * (u_mitt ** 2)) / 2 #designtask2
    delta_T = delta_h / cp #temperatur ratio
    T_t2i = delta_T + T_t2 #Stagnations temperatur efter fläkten
    temp_ratio = T_t2i / T_t2
    FPR = temp_ratio ** (gamma * n / (gamma - 1))# 1.17 med verkningsgrad, 1.32
    P_t2i = P_t2 * FPR
    return [P_t2i, T_t2i, FPR]





