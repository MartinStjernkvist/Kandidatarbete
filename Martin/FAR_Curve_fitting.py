#FAR Curve fitting
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import numpy as np

T = np.linspace(300, 1500, 100)
FAR = [0.034,0.0323,0.031,0.0295,0.0282,0.0267,0.025,0.0227,0.0212,0.020]
T_03 = [500,550,600,650,700,750,820,900,950,1000]

def linjär(x,k,m):
    return k*x+m

fitted_variables = curve_fit(linjär,T_03,FAR)
print(fitted_variables)
plt.plot(T,linjär(T,fitted_variables[0][0],fitted_variables[0][1]), label = 'Fitted curve')
plt.plot(T_03, FAR, label = 'Approximate datapoints from appendix C')
plt.legend()
plt.title('FAR as a function of T_03 for an inlet temperature of 1650[K]')
plt.show()


