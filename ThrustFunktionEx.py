#En funktion som ger Thrust för motorn i exempeluppgiften
import numpy as np

def thrust(H=10668, M0=0.78, OPR=45, FPR=1.48, BPR=12.5, HPCPR=10, TtInlet=1650, fanEpol=0.9, IPCEpol=0.9, HPCEpol=0.9, HPTEpol=0.85, LPTEpol=0.9, CompPloss = 0.04):
    #U på slutet betyder samma som angivet i uppgiften
    HU=10668 #[m]
    M0U=0.78
    OPRU=45
    FPRU=1.48
    BPRU=12.5
    HPCPRU=10
    TtInletU=1650 #[K]
    fanEpolU=0.9
    IPCEpolU=0.9
    HPCEpolU=0.9
    HPTEpolU=0.85
    LPTEpolU=0.9
    CompPlossU = 0.04
    TotalAirMassFlowU = 185 # Detta borde variera med hastighet och höjd
    FARU=0.025 # Fuel air ratio
    R = 287 #[J/kg] gas konstant för luft
    T_0U = T_sls-6.5*10**(-3)*HU
    v0U = M0U*np.sqrt(gamma_a*R*T0U)
    T_sls=288.15 #[K] temperatur för sea level standard
    P_sls = 101325 #[Pa] tryck för sea level standard
    gamma_a = 1.4 #Detta är approximerat som konstant
    g = 9.81 #[m/s^2] Jordens gravitation
    P0U = P_sls * (T_0U/T_sls)^(g*10**3/(6.5*R))
    rho0U = P0U/(R*T_0U) #luftdensitet för uppgiften
    A0=TotalAirMassFlowU/(v0U*rho0U) #[m^2]
    cp_a_sls=10035 #[J/kgK] värmekapacitet för luft sea level standard

    #Saker som ändras med höjd och M0
    IPCPR = OPR/(FPR*HPCPR)
    T_0 = T_sls-6.5*10**(-3)*H #1.29a introkompendium
    P_0 = P_sls * (T_0/T_sls)^(g*10**3/(6.5*R)) #1.29b introkompendium
    a_0 = np.sqrt(gamma_a*R*T_0) #[m/s] ljudhastighet framför motorn
    v_0 = M0*a_0 #[m/s] hastighet på luften framför planet
    rho0 = P0/(R*T_0) #luftdensitet
    TotAirMassFlow = A0*v_0*rho0 #[kg/s] Detta borde vara rätt nu, me  lite osäker
    
    P_02=P_0*(1+((gamma_a-1)/2)*M0**2)**(gamma_a/(gamma_a-1))
    T_02=T_0*(1+((gamma_a-1)/2)*M0**2)

    P_021=P_0*FPR
    T_021=T_02*(P_021/P_0)**((gamma_a-1)/(gamma_a*fanEpol))

    P_026=P_021*IPCPR
    T_026=T_021*(P_026/P_21)**((gamma_a-1)/(gamma_a*IPCEpol))

    P_03=P_026*HPCPR
    T_03=T_026*(P_03/P_26)**((gamma_a-1)/(gamma_a*HPCEpol))

    P_04=P_03*(1-CompPloss)

    FAR = FARU*1 #Fuel air ratio, detta är temp ska fixas senare till en ordentlig sak enligt chart appendix C
    mflow_bypass = TotAirMassFlow/(1+1/BPR)
    mflow_core = TotAirMassFlow - mflow_bypass
    mflow_fuel = mflow_core*FAR

    gamma_g = gamma_a # Detta är temporärt, måste fixas så den är korrekt.
    cp_a = gamma_a*R/(gamma_a-1) # Obs! Detta gäller för ideala gaser
    # cp_f = 0 # specifik värmekapacitet för bränslet
    cp_g = gamma_g*R/(gamma_g-1 )# Obs! Detta gäller för ideala gaser


    T_045 = T_04 - (cp_a*(T_03-T_026))/(cp_g*(1+FAR))#HPC is driven by the HPT
    P_045 = P_04*(T_045/T_04)**((gamma_g)/((gamma_g-1)*HPTEpol))

    Cp_a = cp_a # Behövs dessa???
    Cp_g = cp_g
    # För ideala gaser gäller det att Cp=gamma*n*R/(gamma-1) där n är molsubstans. På så vis får vi att cp=Cp/n
    # I sin tur får vi att Cp_a/Cp_g = cp_a/cp_g! Kanske endast själva förhålllandet som är av intresse?


    T_05 = T_045 - CP_a/(Cp_g*(1+FAR))*((1+BPR)*(T_021-T_02)+(T_026-T_021))# fläkt och IPC drivs av LPT
    P_05 = P_045*(T_05/T_045)**((gamma_g)/((gamma_g-1)*LPTEpol))

    cpr_cold = ((gamma_a+1)/2)**(gamma_a/(gamma_a-1)) #cold critical preasure ratio exit, bypass
    cpr_hot = ((gamma_g+1)/2)**(gamma_g/(gamma_g-1)) #hot critical preasure ratio exit, core

    bypass_choked = False
    core_choked = False

    # beräknas med introkompendie sida 16 vilket ger P_08 aprox = P_05 och P_8 = P_0
    c_8 = np.sqrt(((P_05/P_0)**((gamma_g-1)/gamma_g)-1)*(2/(gamma_g-1))) # TEMPORTÄR MÅSTE BERÄKNAS IFALL CORE EJ CHOKED
    c_18 = np.sqrt(((P_021/P_0)**((gamma_a-1)/gamma_a)-1)*(2/(gamma_a-1))) # TEMPORTÄR MÅSTE BERÄKNAS IFALL CORE EJ CHOKED
    bypass_Thrust = 0
    core_Thrust = 0
    if P_021/P_0 > cpr_cold:
        bypass_choked = True
        P_18 = P_021/cpr_cold # P_021=P_018 aprox
        T_18 = 2/(1+gamma_a)*T_021 
        a_18 = np.sqrt(gamma_a*R*T_18)
        c_18 = a_18 # M18=1 för bypass är choked
        rho_18 = P_18/(R*T_18)
        A_18 = mflow_bypass/(c_18*rho_18)
        bypass_Thrust = A_18*(P_18-P_0)
    if P_05/P_0 > cpr_hot:
        core_choked = True
        P_8 = P_05/cpr_cold # P_08=P_05 aprox
        T_8 = 2/(1+gamma_a)*T_05 
        a_8 = np.sqrt(gamma_g*R*T_8)
        c_8 = a_08 # M18=1 för bypass är choked
        rho_8 = P_8/(R*T_8)
        A_8 = (mflow_core+mflow_fuel)/(c_8*rho_8)
        core_Thrust = A_8*(P_8-P_0)
    F = (mflow_core+mflow_fuel)*c_8+mflow_bypass*c_18-TotAirMassFlow*v_0+core_Thrust+bypass_Thrust #c_8 och c_18 måste beräknas ifall ej choked
    return F