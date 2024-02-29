#%% Symphony 
#-----------------------------------Almänna värden-------------------------------------------------------
import numpy as np

#--------------------------------Standard Värden och Givna Värden------------------------------
T_SLS = 288.15 #kelvin Sea level standard
P_SLS = 101325 #Pa Sea level standard
g = 9.81 #Jordens gravitation
R = 287.05 # luftens unika R
gamma_a = 1.4 # antag att gamma för luft är 1.4
cp_a = gamma_a * R / (gamma_a - 1)  # Obs! Detta gäller för ideala gaser
P_0_60k = 7176 # Pa, linjär interpolering från tabell i introkompendie
M_0 = 1.7 #Detta är cruise punkt
M_2 = 0.6 #Detta var givet enligt en källa

#-------------------------Steg-laster och Relativa Machtal----------------------------------------
steglast_fläkt = 0.6 #** %från handledning %%%%%%%%%%HITTA VÄRDE%%%%%%%%%%%% Detta är stage-load
steglast_HPC = 0.8409 #-5.716+0.00323*EIS --------- EIS = 2030
steglast_LPC = 0.9323 #-8.968 + 0.004877*EIS --------- EIS = 2030
M_rel_fläkt = 1.5 #** %från handledning %%%%%%%%%%%%%%HITTA VÄRDE%%%%%%%%%%%%%
M_rel_HPC = 1.3 #designtask2
M_25 = 0.482 #designtask2 Denna kanska ska förkastas och beräknas utifrån MFP eller liknande?

#-----------------------Effektiviteter--------------------------------
n_fläkt = 0.98 # effektivitet för fläkten
n_inlet = 0.98 # Detta beskriver tryckförlust från P_t1 till P_t2
n_LPC = 0.901 # effektivitet lpc
n_HPC = 0.941 

#-----------------------Approximerat från bild---------------------------
r_tipp_fläkt = 0.915 #approximerat, [m]
r_hub_fläkt = 0.29625 #approximerat, [m]
A_2 = 2.296 # Beräknat i parameter värden dokumentet
r_tipp_LPC = 0.732 # approximerat [m]
r_tipp_HPC = 0.366 #0.4575 approximerat [m]
v_LPC_ingång = 0.630 #hub-tip ratio för LPC
#v_HPC_entry = 0.5613/(1.487-(exp(-0.04286*((FPR*LPC_ratio)+0.5718))))# %%%%%%%%%%%FRÅGA%%%%%%%%%
#v_HPC_entry = 0.85

#------------------Andra antagna värden som behövs optimeras------------------
BPR = 3.5 # förmodligen mellan 3 och 4 ty medium bypass


#%%-------------------------------------Beräkningar---------------------------
T_0, P_0, a_0 = Ambient_air(18288) #Beräknar luftegenskaper för omgivningen
P_0 = P_0_60k

T_t2, P_t2, T_2, P_2 , mdot_luft, a_2, MFP_2 = Station_2(T_0, P_0) #Beräknar massflöde och luftegenskaper innan fläkten

T_t2i, P_t2i, FPR, omega_fläkt = Station_2i(T_t2, P_t2, a_2) #Beräknar luftegenskaper efter fläkten

mdot_bypass = mdot_luft * BPR / (1 + BPR) #bypass flöde
mdot_core_luft = mdot_luft - mdot_bypass #core flöde med endast luft

P_t25, T_t25, LPC_PR = Station_25(P_t2i, T_t2i, omega_fläkt)

P_t3, T_t3 = Station_3(P_t25, T_t25)
#%%---------------------------------Funktioner--------------------------------------
def Ambient_air(h, T_std = T_SLS, P_std = P_SLS, gamma = gamma_a, R = R):
    """
    Ger tillbaka T_0 och P_0 från T och P vid havsnivå
    """
    if h > 11000: # formlerna fungerar ej över 11000m
        T_0 = T_std -6.5*11000/(10**3)
        P_0 = None
    else:  
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
    P_2 = P_t2 / temp2 #osäker ifall detta är rätt
    a_2 = np.sqrt(gamma * R * T_2)
    return [T_t2, P_t2, T_2, P_2, mdot, a_2, MFP]

def Station_2i(T_t2, Pt_2, a_2, cp = cp_a, M_2 = M_2, gamma = gamma_a, M_rel = M_rel_fläkt, r_tipp = r_tipp_fläkt, r_hub = r_hub_fläkt, steglast = steglast_fläkt, n = n_fläkt): # Efter fläkten 
    c_2 = M_2 * a_2 
    u_tipp = np.sqrt((M_rel*a_2) ** 2 - (c_2 ** 2)) #bladhastighet vid tippen
    omega = u_tipp / r_tipp; #vinkelhastighet för fläkt
    rpm = ((omega * 60) / (2 * np.pi)) #rpm för fläkt
    r_mitt = (r_tipp + r_hub) / 2 #radie till berkäning av hastighet u_mitt
    u_mitt = r_mitt * omega #medel av hastighet utöver fläkten 
    delta_h = (steglast * (u_mitt ** 2)) / 2 #designtask2
    delta_T = delta_h / cp #temperatur ratio
    T_t2i = delta_T + T_t2 #Stagnations temperatur efter fläkten
    temp_ratio = T_t2i / T_t2
    FPR = temp_ratio ** (gamma * n / (gamma - 1))# 1.17 med verkningsgrad, 1.32
    P_t2i = P_t2 * FPR
    return [P_t2i, T_t2i, FPR, omega]

def Station_25(P_t2i, T_T2i, omega, v = v_LPC_ingång, r_tipp = r_tipp_LPC, steglast = steglast_LPC, cp = cp_a, n = n_LPC, gamma = gamma_a):
    """
    Ger utloppet av lågtryckskompressorn och/eller inloppet till högtryckskompressorn
    omega = omega_fläkt ty samma axel
    """
    r_hub = v * r_tipp
    r_mitt = (r_tipp + r_hub)/2
    u_mitt = r_mitt * omega
    delta_h = (steglast*((u_mitt**2)*3))/2 #designtask2, gånger 3 pga three stage LPC
    delta_T = delta_h/cp #temperatur ratio
    T_t25 = delta_T + T_t2i #Stagnations temperatur efter fläkten
    temp_ratio = T_t25 / T_t2i
    LPC_PR = temp_ratio ** (gamma * n / (gamma - 1))# 1.17 med verkningsgrad, 1.32
    P_t25 = P_t2i * LPC_PR
    #Om vi har en A_25 kan vi använda MFP för att få M_25 och därmed T_25 och liknande
    return [P_t25, T_t25, LPC_PR]

def Station_3(P_t25, T_t25, M_25 = M_25, M_rel = M_rel_HPC, r_tipp = r_tipp_HPC, steglast = steglast_HPC, cp = cp_a, n = n_HPC, gamma = gamma_a, R = R):
    """
    Utloppet av HPT
    """
    T_25 = T_t25 / (1 + (gamma-1)/2 * M_25 ** 2)
    a_25 = np.sqrt(gamma * R * T_t25)
    c_25 = M_25 * a_25
    u_tipp = np.sqrt((M_rel*a_25) ** 2 - c_25 ** 2)
    omega = u_tipp / r_tipp
    rpm = omega * 60 / (2 * np.pi)
    r_hub = v*r_tipp
    r_mitt = (r_tipp + r_hub)/2
    u_mitt = omega * r_tipp

    delta_h = (steglast*((u_mitt**2)*6))/2 #designtask2, gånger 6 pga three stage LPC
    delta_T = delta_h/cp #temperatur ratio
    T_t3 = delta_T + T_t25 #Stagnations temperatur efter fläkten
    temp_ratio = T_t30 / T_t25
    HPC_PR = temp_ratio ** (gamma * n / (gamma - 1))# 1.17 med verkningsgrad, 1.32
    P_t3 = P_t25 * HPC_PR
    return [P_t3, T_t3]
# %%
