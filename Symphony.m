%% Main_Symphony UTKAST1

clear all
close all
clc
%options = optimset('Display','off');
%------------------FEL VI VET ÄR FEL/KANSKE ÄR FEL------------------------%
% a hastigheterna kanske är fel ty de använda gamma baserat på
% stagnationstemperatur och ej temperatur Detta påverkar Mach_MFP
% funktionen

%Så att alla mixade c_p värden använder f och BPR/masscore på rätt sätt
%samt att dessa löses som optimeringsproblem

% R_mixed kan förändras då temperaturen ändras och har ej korrekta värden
% just nu

% Exhaust_nozzlen är ej klar

% FIXA ENTALPI MIX

%beräkningen av c_p istället för ett "genomsnitt" så blir det värdet i
%början som används, särskillt i förbrännaren är detta viktigt.

%c_p_xm kanske saknar FAR

%stage loading för fläkt? ska det vara 1.8?
%stage loading för turbiner också

%Detta är en viktig förändring!!!

%entalpi i förbränning, med bra referens temperatur, ändra på värden, och implementera sista delen bra 



%Punkter översättningen (kommer ändras)
%2 - innan fläkt
%3 - efter fläkt/innan LPC
%4 - efter LPC/innan HPC
%5 - efter HPC/innan förbränning
%6 - efter förbränning/innan HPT


%c_p: Värmekapacitet vid konstant tryck
%c_v: värmekapacitet vid konstant volym.
%% Input 

a_0 = 295; %speed of sound
c_p = 1005.7; 
T_0 = -56.6+273.15; %kelvin at 60,000ft
gamma_air = 1.4; 
P_0 = 7171.64;
M_0 = 1.7;
c_0 = a_0 * M_0;
EIS = 2030; %entry into service
R = 287.05; %[J/kgK]

%%
%tspace = linspace(0, 1700, 1700);
%plot(tspace, c_p_air(tspace))
%%
%gammalist = zeros(1700);
%for i = 1:1700
%    gammalist(i) = gamma(c_p_air(tspace(i)), R);
%end
%plot(tspace, gammalist)
%ylim([1.25 1.45])
%-----------------------------FAN INPUT VALUES----------------------------%

M_ax_fan = 0.603; %Sida 6 designtask2 Detta borde ej användas
psi_fan = 1.7;%** %Mellan 0.43 och 0.65 sida 10 designtask2
M_rel_fan = 1.5; %från handledning, HITTA NYTT VÄRDE
n_fan = 0.89; %Level 4 tech table 4.4 s. 107 Mattingly 


%Entry, approximerat [m]
r_tip_fan(1) = 0.9525; 
r_hub_fan(1) = 0.309;

%Exit, approximerat [m]
r_tip_fan(2) = 0.9525; 
r_hub_fan(2) = 0.516;

%För att jämföra med uppmätta värden
nu_fan = 44.29/(98.94+exp((0.01850*EIS)-33.31)); %hub-tip ratio för fläkt

%----------THREE STAGE LOW PRESSURE COMPRESSOR (LPC) INPUT VALUES---------%


n_LPC = 0.901; %**Hitt värde för trestegs kompressor
psi_LPC = 0.85;%-8.968 + 0.004877*EIS; %sida 11 designtask2


%Stage 1, approximerat [m]
r_tip_LPC(1) = 0.716;
r_hub_LPC(1) = 0.545;

%Stage 2, approximerat [m]
r_tip_LPC(2) = 0.716;
r_hub_LPC(2) = 0.545;

%Stage 3, approximerat [m]
r_tip_LPC(3) = 0.716;
r_hub_LPC(3) = 0.545;


%För att jämföra med uppmätta värden
nu_LPC_entry = 0.630; %hub-tip ratio steg 1 för LPC
nu_LPC_exit = 0.819; %hub-tip ratio steg 3 för LPC


%----------------------------------MFP------------------------------------%
M_2 = 0.60; %0.6-0.8
A_2 = 2.296; %[m2] estimerat från bild

%--------------HIGH PRESSURE COMPRESSOR (HPC) INPUT VALUES----------------%

M_rel_HPC = 1.3; %designtask2 1.3 från design task 2
M_ax_HPC = 0.482; %designtask2 sida 13 Detta borde en användas
%r_tip_HPC = 0.444; %0.4575; %approximerat [m]
psi_HPC = 0.75;%-5.716+0.00323*EIS; %sida 13 designtask
n_HPC = 0.941;

%Stage 1, approximerat [m]
r_tip_HPC(1) = 0.444;
r_hub_HPC(1) = 0.264; %0.273;

%Stage 2, approximerat [m]
r_tip_HPC(2) = 0.444;
r_hub_HPC(2) = 0.324; %0.334;

%Stage 3, approximerat [m]
r_tip_HPC(3) = 0.444;
r_hub_HPC(3) = 0.349; %0.364;

%Stage 4, approximerat [m]
r_tip_HPC(4) = 0.444;
r_hub_HPC(4) = 0.360; % 0.371;

%Stage 5, approximerat [m]
r_tip_HPC(5) = 0.444;
r_hub_HPC(5) = 0.371; % 0.400;

%Stage 6, approximerat [m]
r_tip_HPC(6) = 0.444;
r_hub_HPC(6) = 0.382; %0.436;



%------------------------COMBUSTION INPUT VALUES--------------------------%
T_t4 = 1800; %Kelvin
n_Combustor = 0.86;
gamma_gas = 1.333; %Introduktionskompendium sida 16
T_tref = 298;

Combustor_ratio = 0.96; %4 procent
Cooling_pre_HPT = 0.1;
Cooling_post_HPT = 0.1;
Combust_massflow_ratio = 1 - Cooling_post_HPT - Cooling_pre_HPT;

%----------SINGLE STAGE HIGH PRESSURE TURBINE (HPT) INPUT VALUES----------%
r_tip_HPT = 0.509;
%r_hub_HPT = 0.444; %%%%OBS ÄNDRA DESSA SEN%%%%%%%
psi_HPT = 3.000;
n_HPT = 0.89; %%%%OBS ÄNDRA DESSA SEN%%%%%%%


%-----------THREE STAGE LOW-PRESSURE TURBINE (LPT) INPUT VALUES-----------%

psi_LPT = -39.26;
n_LPT = 0.90; %%%OBS ÄNDRA DESSA SEN%%%%%

%Stage 1, approximerat [m]
r_tip_LPT(1) = 0.629;

%Stage 2, approximerat [m]
r_tip_LPT(2) = 0.768;

%Stage 3, approximerat [m]
r_tip_LPT(3) = 0.814;

%--------------------------LOBED MIXER INPUT VALUES-----------------------%
A_16 = Area_outer(0) - 1.096;%0.9367; %[m^2] manuellt beräknat från bild
translation = 1; %[m] förskjutning åt höger av spiken
A_6 = 1.096-Area_inner(0,translation); %[m^2] manuellt beräknat från bild, Ska också kunna variera med att spiken flyttas ut och in
n_bypass = 0.98; %Term som beskriver friktionsförluster
A_6A = A_6 + A_16;
%---------------------EXHAUST NOZZLE INPUTVALUES--------------------------%
n_nozzle = 0.98; %Beräkna senare med hjälp av integration?
Exhaust_length = 3.763;
A_8 = Area_exhaust(Exhaust_length,translation);%1.18; %meter beräknat från bild

%% Pre calculations (inlet)

P_t2 = P_0*(((1+((gamma_air-1)/2)*M_0^2))^(gamma_air/(gamma_air-1))); %Mattingly

T_t2 = T_0*((1+(((gamma_air-1)/2)*M_0^2)));
T_2 = T_t2/((1+(((gamma_air-1)/2)*M_2^2)));

Entalphy_above_T_tref = (c_p_air(T_t2)+c_p_air(T_tref))/2*(T_t2-T_tref);

c_p_2 = c_p_air(T_t2); % Stagnations temperatur används för att man kan tänka sig att gasen bromsas till c = 0, värms upp och sedan höjs hastigheten igen
gamma_2 = gamma(c_p_2,R);

%% Mass Flow Parameter (MFP)

MFP = M_2*(sqrt(gamma_2/R))*(1+(((gamma_2-1)/2)*M_2^2))^((gamma_2+1)/(2*(1-gamma_2)));

massflow = (MFP*P_t2*A_2)/(sqrt(T_t2)); %Mattingly (1.3)


%-------------------------BYPASS AND MASSFLOW-----------------------------%
%r_ratio_BPR = r_tip_fan(1)/r_tip_LPT(3);
%BPR = ((log(1.9176-(r_ratio_BPR*1.25)))/(-0.2503))-0.6410;
BPR = 4.5;
massflow_bypass = massflow * BPR/(1+BPR);
massflow_core =  massflow - massflow_bypass;

%% FPR  

%Jämförelse av hubben genom hub tip ratio och mätta värden.
r_hub_fan_nu = nu_fan*r_tip_fan(1);
ratio_hub_fan = (1-(r_hub_fan_nu/r_hub_fan(1)))*100; %Procentuell skillnad

%Axial hastighet + speed of sound innan fläkt
a2 = sqrt(gamma_2*R*T_2); %speed of sound

M_ax_fan = M_2; %-------------------------------------Kom ihåg denna ändringen
Ca_fan = M_ax_fan*a2; 

%Bladhastigheter
u_fan = sqrt((M_rel_fan*a2)^2-(Ca_fan^2)); %bladhastighet vid tipp

omega_fan = u_fan/r_tip_fan(1); %vinkelhastighet för fläkt

rpm_fan = ((omega_fan*60)/(2*pi)); %rpm för fläkt

%Tar inte hänsyn till konisk form på hubben
r_mid_fan = (r_tip_fan(1) + r_hub_fan(1))/2; %radie till berkäning av hastighet u_mid

u_mid_fan = r_mid_fan*omega_fan; %medel av hastighet utöver fläkten 


%Steglast och flödesfaktor
delta_h_fan = (psi_fan*(u_mid_fan^2))/2; %designtask2
delta_T_fan = Temp_change(T_t2, delta_h_fan, 1, 0);%delta_T_fan = delta_h_fan/c_p_2; %temperatur ratio


%Beräkning av FPR
T_t21 = delta_T_fan + T_t2; %temperaturen efter fläkten

temp_ratio_fan = T_t21/T_t2;

FPR = (temp_ratio_fan)^(gamma_2*n_fan/(gamma_2-1)); %1.17 med verkningsgrad, 1.32

P_t21 = P_t2*FPR;

c_p_21 = c_p_air(T_t21);
gamma_21 = gamma(c_p_21,R);

A_21 = pi*r_tip_LPC(1)^2*(1-(r_hub_LPC(1)/r_tip_LPC(1))^2);

MFP_21 = massflow_core*sqrt(T_t21)/(P_t21*A_21);
M_21 = fsolve(@(M) Mach_MFP(M, gamma_21, R) - MFP_21, 0);

%%
%Mspace = linspace(0,2,1000);
%plot(Mspace,Mach_MFP(Mspace,gamma_21,R)-MFP_21)

%% Three Stage Low pressure compressor (LPC) 

%samma axel som fläkten
omega_LPC = omega_fan; %pga samma axel, alltså samma vinkelhastighet.

%u_mid efter mått i bild
r_mid_LPC = (r_tip_LPC + r_hub_LPC)./2; %radie till berkäning av hastighet u_mid

u_mid_LPC = r_mid_LPC.*omega_LPC; %medel av hastighet utöver fläkten 

%u_mid linjär approximation - bara tre steg.

%steglast och flödesfaktor

u_sum_LPC = (u_mid_LPC(1))^2+(u_mid_LPC(2))^2+(u_mid_LPC(3))^2; %Three stage LPC

delta_h_LPC = (psi_LPC*((u_sum_LPC)))/2; %designtask2
delta_T_LPC = Temp_change(T_t21, delta_h_LPC, 1, 0);%delta_T_LPC = delta_h_LPC/c_p_21; %temperatur ratio

%beräkning av LPC ratio
T_t25 = delta_T_LPC + T_t21; %temperaturen efter LPC

temp_ratio_LPC = T_t25/T_t21;

LPC_ratio = (temp_ratio_LPC)^(gamma_21*n_LPC/(gamma_21-1)); %1.17 med verkningsgrad, 1.32

P_t25 = P_t21*LPC_ratio;

c_p_25 = c_p_air(T_t25);
gamma_25 = gamma(c_p_25,R);

A_25 = r_tip_HPC(1)^2*pi*(1-(r_hub_HPC(1)/r_tip_HPC(1))^2); % Detta är helt fel eller något xdddd

MFP_25 = massflow_core*sqrt(T_t25)/(P_t25*A_25);
M_25 = fsolve(@(M) Mach_MFP(M, gamma_25, R) - MFP_25, 0);

T_25 = T_t25/((1+(((gamma_25-1)/2)*M_25^2)));

%% Six stage High Pressure Compressor (HPC)

%Axial hastighet + speed of sound innan fläk
a25 = sqrt(gamma_25*R*T_25); %--------------
M_ax_HPC = M_25; % --------------------------------------- Kom ihåg detta
Ca_HPC = M_ax_HPC*a25; 

%Bladhastigheter
u_HPC = sqrt((M_rel_HPC*a25)^2-(Ca_HPC^2)); %bladhastighet vid tipp

omega_HPC = u_HPC/r_tip_HPC(1); %vinkelhastighet för fläkt

rpm_HPC = ((omega_HPC*60)/(2*pi)); %rpm för fläkt

%u_mid efter mått i bild
r_mid_HPC = (r_tip_HPC + r_hub_HPC)./2; %radie till berkäning av hastighet u_mid

u_mid_HPC = r_mid_HPC.*omega_HPC; %medel av hastighet utöver fläkten 

u_sum_HPC = sum(u_mid_HPC.^2);

%steglast och flödesfaktor
delta_h_HPC = (psi_HPC*(u_sum_HPC))/2; %designtask2, 6 stage???
delta_T_HPC = Temp_change(T_t25, delta_h_HPC, 1, 0); %delta_T_HPC = delta_h_HPC/c_p_25; %temperatur ratio

%beräkning av HPC ratio
T_t3 = delta_T_HPC + T_t25; %temperaturen efter HPC

temp_ratio_HPC = T_t3/T_t25;

HPC_ratio = (temp_ratio_HPC)^(gamma_25*n_HPC/(gamma_25-1)); %1.17 med verkningsgrad, 1.32

P_t3 = P_t25*HPC_ratio;
c_p_3 = c_p_air(T_t3);
gamma_3 = gamma(c_p_3,R);

A_3 = r_tip_HPC(6)^2*pi*(1-(r_hub_HPC(6)/r_tip_HPC(6))^2); % Detta är helt fel eller något xdddd

MFP_3 = massflow_core*sqrt(T_t3)/(P_t3*A_3);
M_3 = fsolve(@(M) Mach_MFP(M, gamma_3, R) - MFP_3, 0);

% Tillägg ev. kylning. ca 20% Ta enbart ut kylning efter HPC
%sida 104-105 Mattingly.

%% Combustor 

%NYTT GAMMA OCH CP FÖR GAS OCH AIR
R_mixed = 380; %Fixa detta ordentligt ööööööööööööööööööööööö---------------------------------------------------
LHV = 43.5*10^6; %Heat value: Joule/kg för bränslet, detta är lite över den undre gränsen för SAF, borde vara ett rimligt värde, kan behöva dubbelkollas
FAR_guess = 0.02; % 1/(1+BPR) kommer att minska ty endas:              - (T_t3)*c_p_3*Combust_massflow_ratio
FAR = fsolve(@(f) Entalpi_mix(f, T_t4, T_tref, LHV, Combust_massflow_ratio) - (delta_h_HPC - delta_h_LPC - delta_h_fan - Entalphy_above_T_tref)*Combust_massflow_ratio, FAR_guess); %borde vara gåner 10^3 ty cp är per g och ej per kg
%% Allting är helt fuckat, jag hatar livet
%flist = linspace(0,0.1,2000);
%entalpilist = zeros(2000);
%for i = 1:2000
%    entalpilist(i) = Entalpi_mix(flist(i), T_t4, T_0, LHV, Combust_massflow_ratio) - delta_h_HPC - delta_h_LPC - delta_h_fan;
%end
%plot(flist, entalpilist)
%%
c_p_4 = c_p_mixed(T_t4,Combust_massflow_ratio, FAR);
%entalpi
massflow_hot = massflow_core*(1+FAR);

gamma_4 = gamma(c_p_4, R_mixed);
P_t4 = P_t3*Combustor_ratio; %KOLLA SÅ RÄTT RATIO, tror rätt att trycket minskar?

%% (High Pressure Turbine) (HPT)

omega_HPT = omega_HPC; %samma vinkelhastighet

%CHECK AN? BLADE ROOT STRESS LEVELS.

%FINNS INGET NU_HPT???????????????????????
%%%%%%r_hub_HPT = v_HPT_entry*r_tip_HPT; %beräkning av r_hub_HPT

%SKA INTE VARA U MID UTAN BARA U, SE SIDA 16
%r_mid_HPT = (r_tip_HPT + r_hub_HPT)/2; %radie till berkäning av hastighet u_mid

% mixning som är lite fel antar jag
%c_p_4m = ( c_p_4*Combust_massflow_ratio+c_p_3*Cooling_pre_HPT ) / (Combust_massflow_ratio + Cooling_pre_HPT);
%T_t4m = (c_p_4*T_t4 * (1/(1+BPR)*Combust_massflow_ratio +FAR) + c_p_3*T_t3*(1/(1+BPR))*Cooling_pre_HPT ) / ((1/(1+BPR)*(1-Cooling_post_HPT) + FAR)*c_p_4m);
T_t4m = fsolve(@(T_t) Mixer(T_t, T_t4, c_p_4,Combust_massflow_ratio, FAR, T_t3, c_p_3, Cooling_pre_HPT), T_t4); %Detta är efter kylning inblandat
c_p_4m = c_p_mixed(T_t4m, 1-Cooling_post_HPT, FAR);
gamma_4m = gamma(c_p_4m, R_mixed);

u_HPT = r_tip_HPT*omega_LPC; %medel av hastighet utöver fläkten 

%steglast och flödesfaktor
delta_h_HPT = delta_h_HPC / (1-Cooling_post_HPT + FAR);% division se ex uppgift ty iolika massflöden %(psi_HPT*(u_HPT^2)); %designtask2, bara ett steg, annan formel för
%enligt stage loading på sida 16 i designtask 2.
delta_T_HPT = Temp_change(T_t4m, delta_h_HPT, 1-Cooling_post_HPT, FAR);%delta_T_HPT = delta_h_HPT/c_p_4m; %temperatur ratio

%beräkning av HPT
T_t45 = T_t4 - delta_T_HPT; %temperaturen över HPT

temp_ratio_HPT = T_t4/T_t45;

HPT_ratio = (temp_ratio_HPT)^(gamma_4m/((gamma_4m-1)*n_HPT)); %1.17 med verkningsgrad, 1.35

P_t45 = P_t4/HPT_ratio;
c_p_45 = c_p_mixed(T_t45,1-Cooling_post_HPT, FAR);
gamma_45 = gamma(c_p_45, R_mixed);

%c_p_45m = c_p_45*(1-Cooling_post_HPT)+c_p_3*Cooling_post_HPT; %fel, inkluderar ej far, se andra mixed cp
%T_t45m = (c_p_45*T_t45 * (1/(1+BPR)*(1-Cooling_post_HPT) +FAR) + c_p_3*T_t3*(1/(1+BPR))*Cooling_post_HPT ) / ((1/(1+BPR) + FAR)*c_p_45m);
T_t45m = fsolve(@(T_t) Mixer(T_t, T_t45, c_p_45,1-Cooling_post_HPT, FAR, T_t3, c_p_3, Cooling_post_HPT), T_t45); %Detta är efter kylning inblandat
c_p_45m = c_p_mixed(T_t45m, 1, FAR);
gamma_45m = gamma(c_p_45m, R_mixed);

%% THREE STAGE LOW-PRESSURE TURBINE (LPT)

%samma axel som fläkten
omega_LPT = omega_fan; %pga samma axel, alltså samma vinkelhastighet.

%u_mid efter mått i bild
%r_mid_LPC = (r_tip_LPC + r_hub_LPC)./2; %radie till berkäning av hastighet u_mid

u_LPT = r_tip_LPT.*omega_LPT; %medel av hastighet utöver fläkten 
u_sum_LPT = sum(u_LPT.^2);

%steglast och flödesfaktor

delta_h_LPT = (delta_h_LPC + delta_h_fan*(1+BPR))/(1+FAR); %(psi_LPC*(u_sum_LPT)); %designtask2, gånger 3 pga three stage LPC
delta_T_LPT = Temp_change(T_t45m, delta_h_LPT, 1, FAR);%delta_T_LPT = delta_h_LPT/c_p_45m; %temperatur ratio
% Samma för delta_h_LPT som för HPT kommentar
% ------------------------------------------------------------------------------------------------------------------------------
%beräkning av LPC ratio
%beräkning av HPT
T_t5 = T_t45m - delta_T_LPT; %temperaturen över HPT

temp_ratio_LPT = T_t45m/T_t5;

LPT_ratio = (temp_ratio_LPT)^(gamma_45m/((gamma_45m-1)*n_LPT)); %1.17 med verkningsgrad, 1.35

P_t5 = P_t45/LPT_ratio;

c_p_5 = c_p_mixed(T_t5, 1, FAR);
gamma_5 = gamma(c_p_5, R_mixed);

%% Lobbade mixern
%Mspace = linspace(0, 3, 100);
%plot(Mspace,real(Mach_MFP(Mspace, gamma_5, R_mixed)) - MFP_6)
%%
P_t16 = P_t21*n_bypass;
MFP_16 = massflow_bypass*sqrt(T_t21)/(P_t16*A_16); %T_t21 och P_t21 används ty antas "isentropiskt".
M_16 = fsolve(@(M) Mach_MFP(M, gamma_21, R) - MFP_16, 0);
MFP_6 = massflow_hot*sqrt(T_t5)/(P_t5*A_6); %T_t5 och P_t5 används ty antags isentropisk
M_6 = fsolve(@(M) real(Mach_MFP(M, gamma_5, R_mixed)) - MFP_6, 0.5);

massflow_exhaust = massflow_hot + massflow_bypass;

% mixern är lobbad 
R_6A = 330; % mindre bränsle så mer likt luft? Denna ska beräknas analytiskt
%c_p_6A = (c_p_5*massflow_hot + c_p_21*massflow_bypass) / (massflow_hot+massflow_bypass); %------Dessa kanske löses som ett optimeringsproblem????
%T_t6A = (c_p_21 * T_t21 * massflow_bypass + c_p_5 * T_t5 * massflow_hot) / ((massflow_exhaust)*c_p_6A);
T_t6A = fsolve(@(T_t) Mixer(T_t, T_t5, c_p_5, 1, FAR, T_t21, c_p_21, BPR), (T_t5+T_t21)/2);
c_p_6A = c_p_mixed(T_t6A,1+BPR,FAR);
gamma_6A = gamma(c_p_6A, R_6A);

P_6 = P_t5 / ((1+((gamma_5-1)/2)*M_6^2)^(gamma_5/(gamma_5-1)));
P_16 = P_t16 / ((1+((gamma_21-1)/2)*M_16^2)^(gamma_21/(gamma_21-1)));

I_in = P_6*A_6*(1+gamma_5*M_6^2) + P_16*A_16*(1+gamma_21*M_16^2);
P_t6A = (P_t5*massflow_hot+P_t16*massflow_bypass) / (massflow_exhaust); % Detta kanske är fel

M_6A = fsolve(@(M) Impuls_ut(M, P_t6A, A_6A, gamma_6A) - I_in, 0.6); % fixa P_t_6A eller någonting!
T_6A = T_t6A/((1+(((gamma_6A-1)/2)*M_6A^2)));
a_6A = sqrt(gamma_6A*R_6A*T_6A);
c_6A = M_6A*a_6A;
%%
%Mspace = linspace(0, 3, 100);
%Impuls = zeros(100);
%for i = 1:100
%    Impuls(i) = Impuls_ut(Mspace(i), P_t6A, A_6A, gamma_6A) - I_in;
%end
%plot(Mspace,Impuls)

%% Exhaust nozzle
%Detta är extremt primitivt och fel!
T_t8 = T_t6A;
P_t8 = P_t6A*n_nozzle;
c_p_8 = c_p_6A;
gamma_8 = gamma_6A;
R_8 = R_6A;
MFP_8 = massflow_exhaust * sqrt(T_t8)/(P_t8*A_8);
M_8 = fsolve(@(M) real(Mach_MFP(M, gamma_8, R_8)) - MFP_8, 2);
T_8 = T_t8 / ((1+(((gamma_8-1)/2)*M_8^2)));
a_8 = sqrt(gamma_8*R_8*T_8);
c_8 = M_8*a_8;
P_8 = P_t8 / ((1+((gamma_8-1)/2)*M_8^2)^(gamma_8/(gamma_8-1)));
Thrust = massflow_exhaust*c_8-massflow*c_0+A_8*(P_8-P_0);
disp(Thrust)
%% alternativ till exhausten
%användet dVdx och ode45
% T_t8 = T_t6A;
% iterations = Exhaust_length*1000+1;
% x_list = linspace(0,Exhaust_length,iterations);
% c_8 = zeros(1,iterations);
% c_8(1) = c_6A;
% M_8 = zeros(1,iterations);
% M_8(1) = M_6A;
% 
% for i = 2:iterations
%     A_8 = Area_exhaust(x_list(i),translation);
%     j=i-1;
%     %V = [x0,c0,M0,Tt,R,g,translation]
%     V = dVdx( x_list(i), x_list(j),[c_8(j),M_8(j), T_t8, R_6A, gamma_6A, translation] );
%     c_8(i) = V(1);
%     M_8(i) = V(2);
% end
% P_8 = P_t8 / ((1+((gamma_8-1)/2)*M_8(end)^2)^(gamma_8/(gamma_8-1)));
% Thrust = massflow_exhaust*c_8(end)-massflow*c_0+A_8*(P_8-P_0);
%%
%plot (x_list,Area_exhaust(x_list,translation))
%%
%plot (x_list,c_8)
%%
%plot (x_list,M_8)
%%
%Mspace = linspace(0,3, 100);
%plot (Mspace,real(Mach_MFP(Mspace, gamma_21, R)) - MFP_16)
%% OUTPUT

OPR = FPR * LPC_ratio * HPC_ratio;
fprintf('BPR = %0.2f\nMassflow = %0.2f kg/s\n',BPR,massflow)

%%
function MFP = Mach_MFP(M, gamma, R)
    MFP = (1+(gamma-1)/2*M.^2).^((gamma+1)/(2*(1-gamma))).*M.*sqrt(gamma/R);
end

function c_p = c_p_air(T)
    A0 = 1047.63;
    A1 = -0.39;
    A2 = 8.89 * 10^-4;
    A3 = -1.64*10^-7;
    A4 = -6.65*10^-10;
    A5 = 6.03*10^-13;
    A6 = -2.07*10^-16;
    A7 = 2.59*10^-20;
    c_p = A0 + A1.*T + A2.*T.^2 + A3.*T.^3 + A4.*T.^4 + A5.*T.^5 + A6.*T.^6 + A7.*T.^7;
end

function c_p = c_p_fuel(T)
    A0 = 309.08;
    A1 = 9.24;
    A2 = -1.87*10^-2;
    A3 = 2.43*10^-5;
    A4 = -1.85*10^-8;
    A5 = 8.08*10^-12;
    A6 = -1.9*10^-15;
    A7 = 1.89*10^-19;
    c_p = A0 + A1.*T + A2.*T.^2 + A3.*T.^3 + A4.*T.^4 + A5.*T.^5 + A6.*T.^6 + A7.*T.^7;
end

function g = gamma(c_p,R)
    g = c_p/(c_p-R); % ------------------------------- Gör om så att R är ett input värde ty R för mix är skild från R för luft
end

function I = Impuls_ut(M, Pt, A, gamma)
    I = Pt*A*(1+gamma*M^2) / ((1+((gamma-1)/2)*M^2)^(gamma/(gamma-1))); %Pt/P och impuls bevaring
end

function entalpi_mix = Entalpi_mix(f, T_t4, T_t3, LHV, Combust_massflow_ratio)
    % c_p_3_mixed = c_p_mixed(T_t3,Combust_massflow_ratio,f);
    % c_p_4 = c_p_mixed(T_t4,Combust_massflow_ratio,f);
    % c_p_g = (c_p_3_mixed+c_p_4)/2; %linjär interpolering för genomsnittlig cp för mixen?
    entalpi_after = 0;
    iterations = 100;
    deltaT = (T_t4-T_t3)/iterations;
    for i = 1:iterations
        entalpi_after = entalpi_after + deltaT*c_p_mixed(T_t3+deltaT*(i-0.5)/100,Combust_massflow_ratio,f);
    end
    entalpi_mix = entalpi_after*(Combust_massflow_ratio + f) - LHV*f; % ----- Borde vara gånger 10^3 ty c_p är per g och ej per kg
end

function c_p = c_p_mixed(T_t, m, mf)
    c_p = (c_p_air(T_t)*m + c_p_fuel(T_t)*mf)/(m+mf);
end

function enthalpy_difference = Mixer(T_t, T_t1, c_p1, m1, f1, T_t2, c_p2, m2)
    %f1 mängden bränsle, så m1+f1 är totalt massflöde för varmt
    m = m1 + m2;
    enthalpy_before = c_p1*T_t1*(m1+f1)+c_p2*T_t2*m2;
    enthalpy_after = c_p_mixed(T_t,m,f1)*T_t*(m+f1);
    enthalpy_difference = enthalpy_after - enthalpy_before;
end

function V1 = dVdx(x1,x0,V)
    %V = [x0,c0,M0,Tt,R,g,translation]
    c0 = V(1);
    M0 = V(2);
    Tt = V(3);
    R = V(4);
    g = V(5);
    t = V(6);
    A1 = Area_exhaust(x1,t);
    A0 = Area_exhaust(x0,t);
    a0 = c0/M0;
    dA = A1-A0;
    dc = dA/A0*c0/(M0^2-1);
    c1 = c0 + dc;
    M1 = c1/a0;
    T1 = Tt / (1+(g-1)/2*M1^2);
    %a1 = sqrt(g*R*T1);
    V1(1) = c1;
    V1(2) = M1;
    %Lös detta ekv systemet, detta ska egentligen göras men öööööö
    %a1 = sqrt(g*R*T1);
    %T1 = Tt / (1+(gamma-1)/2*M1^2);
    %M1 = c1 / a1;
    
    
end

function area = Area_inner(x,translation)
    %translation är förflyttning av center spiken högerut, variabler är i
    %cm antar jag?
    x = x - translation;
    if x <=0
        radie = 0.361;
    elseif x <= 1.04 %temporärt
        radie = 0.361 - 0.122 * x/1.04;
    elseif x <= 1.91
        radie = 0.239;
    elseif x <= 3.76
        radie = 0.239 + 0.087*(x-1.91)/1.85;
    elseif x <= 5.61  % slutlängd är x = 5.61 meter
        radie = 0.326-0.326*(x-3.76)/1.85;
    else
        radie = 0;
    end
    area = pi.*radie.^2;
end

function area = Area_outer(x)
   %x är noll vid något ställe xdd
   %cm antar jag?
    if x <=0
        radie = 0.805;
    elseif x <= 0.652 %temporärt
        radie = 0.805 - 0.153 * x/0.652;
    elseif x <= 2.502
        radie = 0.652;
    elseif x <= 3.263
        radie = 0.652 + 0.065*(x-2.502)/0.759;
    elseif x <= 3.763  % slutlängd är x = 3.763 meter = exhaust length
        radie = 0.717+0.022*(x-3.263)/0.5;
    else 
        radie = NaN;
    end
    area = pi.*radie.^2;
end

function area = Area_exhaust(x,translation)
    area = Area_outer(x) - Area_inner(x, translation);
end

function deltaT = Temp_change(T_t0, delta_H, m, f)
    steps = 500;
    h = delta_H/steps;
    T_t = T_t0;
    for i = 1:steps
        c_p = c_p_mixed(T_t,m,f);
        T_t = h/c_p + T_t;
    end
    deltaT = T_t-T_t0;
end