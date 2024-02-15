#temp
import numpy as np

P0 = 23860 
M0 = 0.78
gamma = 1.4
T0 = 218.81
Tt2 = (1+(gamma-1)/2*M0**2)*T0 #ex uppgift facit
Pt2 = P0*(1+(gamma-1)/2*M0**2)**(gamma/(gamma-1))
g = 9.81
R = 287
MFP0 = M0*np.sqrt(gamma*g/R)*(1+(gamma-1)/2*M0**2)**((gamma+1)/(2*(1-gamma))) # 1.4 mattingly
mdot = 185
A0 = mdot*np.sqrt(Tt2)/(Pt2*MFP0)
print(A0)