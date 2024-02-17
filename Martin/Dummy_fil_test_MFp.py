# temp
import numpy as np
from FORMLER_Mattingly import MFP

"""
Konstanter
"""
P0 = 23860
M0 = 0.78
M2 = 0.42  # mach tal precis innan motorn
gamma = 1.4
T0 = 218.81
R = 287
mdot_ex = 185  # inflöde fån exempeluppgiften
P_effektivitet = 0.99  # Pt2 Verklig / Pt2 teoretisk

"""
Motor dimensioner
"""
diameter = 2.06
radius = diameter / 2
A0_verklig = np.pi * radius ** 2 * (
            1 - (2.2 / 9.4) ** 2)  # estimering baserat på bild och att fläkt diameter är 2.06


Tt2 = (1 + (gamma - 1) / 2 * M0 ** 2) * T0  # ex uppgift facit
Pt2 = P0 * (1 + (gamma - 1) / 2 * M0 ** 2) ** (gamma / (gamma - 1)) * P_effektivitet

MFP0 = M2 * np.sqrt(gamma / R) * (1 + (gamma - 1) / 2 * M2 ** 2) ** ((gamma + 1) / (2 * (1 - gamma)))  # 1.3 mattingly

A0 = mdot_ex * np.sqrt(Tt2) / (Pt2 * MFP0)
mdot_test = A0_verklig * (Pt2 * MFP0) / np.sqrt(Tt2)

print(f'A0: {A0}')
print(f'A0_verklig: {A0_verklig}')
print(f'mdot_test:{mdot_test}')
