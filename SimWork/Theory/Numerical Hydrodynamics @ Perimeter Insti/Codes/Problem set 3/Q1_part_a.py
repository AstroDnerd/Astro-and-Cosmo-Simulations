#function u(x) = sin(x)
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
del_x = np.array([1,0.1,0.01,0.001])
log_del_x = np.zeros(len(del_x))
#fit function
def lin(x,a,b):
    return a*x+b
#forward difference
D_fwd_u = np.zeros(len(del_x))
trunc_error_fwd = np.zeros(len(del_x))
for x in range(len(del_x)):
    D_fwd_u[x] = (np.sin(1+del_x[x])-np.sin(1))/del_x[x]
    trunc_error_fwd[x] = np.log10(np.cos(1)-D_fwd_u[x])
    log_del_x[x] = np.log10(del_x[x])
popt, pcov = curve_fit(lin, log_del_x, trunc_error_fwd)
plt.scatter(log_del_x,trunc_error_fwd)
plt.plot(log_del_x, lin(log_del_x, *popt), 'r-',
         label='fit: a=%5.3f, b=%5.3f' % tuple(popt))
plt.title("Forward Difference")
plt.xlabel("log(delta_x)")
plt.ylabel("log(truncation_error)")
plt.legend()
plt.show()
#centered difference
D_cntr_u = np.zeros(len(del_x))
trunc_error_cntr = np.zeros(len(del_x))
for x in range(len(del_x)):
    D_cntr_u[x] = (np.sin(1+del_x[x])-np.sin(1-del_x[x]))/(2*del_x[x])
    trunc_error_cntr[x] = np.log10(np.cos(1)-D_cntr_u[x])
popt, pcov = curve_fit(lin, log_del_x, trunc_error_cntr)
plt.scatter(log_del_x,trunc_error_cntr)
plt.plot(log_del_x, lin(log_del_x, *popt), 'r-',
         label='fit: a=%5.3f, b=%5.3f' % tuple(popt))
plt.title("Centered Difference")
plt.xlabel("log(delta_x)")
plt.ylabel("log(truncation_error)")
plt.legend()
plt.show()

#Value of a is close to fwd in first and close to 2 in cntr
