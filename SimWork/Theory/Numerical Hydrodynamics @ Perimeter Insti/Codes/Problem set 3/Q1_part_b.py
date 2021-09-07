#function u(x) = sin(x)
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
del_x = np.array([1,0.1,0.01,0.001])
log_del_x = np.zeros(len(del_x))
#fit function
def lin(x,a,b):
    return a*x+b
#backward difference
D2_bck_u = np.zeros(len(del_x))
trunc_error_bck = np.zeros(len(del_x))
for x in range(len(del_x)):
    D2_bck_u[x] = (2*np.sin(1)-5*np.sin(1-del_x[x])+4*np.sin(1-2*del_x[x])-np.sin(1-3*del_x[x]))/(del_x[x]**2)
    trunc_error_bck[x] = np.log10(np.abs(-np.sin(1)-D2_bck_u[x]))
    log_del_x[x] = np.log10(del_x[x])
popt, pcov = curve_fit(lin, log_del_x, trunc_error_bck)
plt.scatter(log_del_x,trunc_error_bck)
plt.plot(log_del_x, lin(log_del_x, *popt), 'r-',
         label='fit: a=%5.3f, b=%5.3f' % tuple(popt))
plt.title("Backward Difference for second derivative")
plt.xlabel("log(delta_x)")
plt.ylabel("log(truncation_error)")
plt.legend()
plt.show()
