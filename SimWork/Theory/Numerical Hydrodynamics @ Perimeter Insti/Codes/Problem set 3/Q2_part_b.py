#heat equation scheme (1st order fwd in t and 2nd ordered cntr in x
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
del_x_array = np.array([0.02,0.01,0.005])
del_t_array = np.array([0.2,0.1,0.005])
nu_d = 0.001
t_start = 1
t_end = 3
x_start = 0
x_end = 1
#analytical solution
def ana_sol(x,t):
    sol = 0
    for n in range(1,51):
        sol+= (2*np.sin(np.pi*n/2)*np.sin(np.pi*n*x)*np.exp(-np.pi*np.pi*n*n*nu_d*t))
    return sol

#numerical scheme
def u_i_nplus1(delt,delx,nu_d,u_i_n,u_iplus1_n,u_iminus1_n):
    return (u_i_n + ((delt*nu_d)/(delx**2))*(u_iplus1_n-2*u_i_n+u_iminus1_n))

#fractional global error
def E_i_n(u_i_n,v_i_n):
    return (np.abs(u_i_n-v_i_n)/(np.abs(v_i_n)))

#Case 1, del_x = 0.002 and del_t varies
del_x = del_x_array[0]
x_size =  1+int((x_end - x_start)/del_x)
x_array = np.linspace(x_start,x_end,x_size)

#initial condition
v_x_1 = np.zeros(x_size)
for x in range(x_size):
    v_x_1[x] = ana_sol(x_array[x],1)

#final analytical result 
v_x_3 = np.zeros(x_size)
for x in range(x_size):
    v_x_3[x] = ana_sol(x_array[x],3)

#evolution
for del_t in del_t_array:
    t_size = int((t_end-t_start)/del_t)
    t_array = np.linspace(t_start+del_t,t_end,t_size)
    evolution_array = np.zeros([t_size+1,x_size])
    evolution_array[0,:] = v_x_1
    for t in range(1,t_size+1):
        for x in range(x_size):
            u_i_n = evolution_array[t-1,x]
            if x == x_size-1:
                u_iplus1_n = evolution_array[t-1,0]
            else:
                u_iplus1_n = evolution_array[t-1,x+1]
            if x ==0:
                u_iminus1_n = evolution_array[t-1,1]
            else:
                u_iminus1_n = evolution_array[t-1,x-1]
            evolution_array[t,x] = u_i_nplus1(del_t,del_x,nu_d,u_i_n,u_iplus1_n,u_iminus1_n)
    u_x_3 = evolution_array[t_size,:]
    #error
    glob_err = np.zeros(x_size)
    for x in range(x_size):
        glob_err[x] = E_i_n(u_x_3[x],v_x_3[x]+10e-5)
    #excluding end points coz analytical solution is 0 there
    plt.scatter(x_array[1:-1],glob_err[1:-1])
    plt.title("(i) Global fraction error for del_x = %5.3f, del_t = %5.3f"% (del_x,del_t))
    plt.show()

    
#Case 2, del_t = 0.2 and del_x varies
del_t = del_t_array[0]
for del_x in del_x_array:
    x_size =  1+int((x_end - x_start)/del_x)
    x_array = np.linspace(x_start,x_end,x_size)

    #initial condition
    v_x_1 = np.zeros(x_size)
    for x in range(x_size):
        v_x_1[x] = ana_sol(x_array[x],1)

    #final analytical result 
    v_x_3 = np.zeros(x_size)
    for x in range(x_size):
        v_x_3[x] = ana_sol(x_array[x],3)

    #evolution

    t_size = int((t_end-t_start)/del_t)
    t_array = np.linspace(t_start+del_t,t_end,t_size)
    evolution_array = np.zeros([t_size+1,x_size])
    evolution_array[0,:] = v_x_1
    for t in range(1,t_size+1):
        for x in range(x_size):
            u_i_n = evolution_array[t-1,x]
            if x == x_size-1:
                u_iplus1_n = evolution_array[t-1,0]
            else:
                u_iplus1_n = evolution_array[t-1,x+1]
            if x ==0:
                u_iminus1_n = evolution_array[t-1,1]
            else:
                u_iminus1_n = evolution_array[t-1,x-1]
            evolution_array[t,x] = u_i_nplus1(del_t,del_x,nu_d,u_i_n,u_iplus1_n,u_iminus1_n)
    u_x_3 = evolution_array[t_size,:]
    #error
    glob_err = np.zeros(x_size)
    for x in range(x_size):
        glob_err[x] = E_i_n(u_x_3[x],v_x_3[x]+10e-5)
    #excluding end points coz analytical solution is 0 there
    plt.scatter(x_array[1:-1],glob_err[1:-1])
    plt.title("(ii) Global fraction error for del_t = %5.3f, del_x = %5.3f"% (del_t,del_x))
    plt.show()

