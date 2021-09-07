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
def u_i_nplus1(x_start,x_end,delt,delx,nu_d,u_n):
    '''writing the CN scheme in matrix form, we have:
        Lu_{n+1} = Qu_{n}
        Where L and Q are matrices of the form:
        L = ((1+lambda)-(lambda/2)*(S_{+} + S_{-})
        Q = ((1-lambda)+(lambda/2)*(S_{+} + S_{-})
        and lambda = (del_t*nu_d)/((del_x)^2)
        S_{+} is of the form:
        |01000.....0|
        |00100.....0|
        |00010.....0|
        |...........|
        |........100|
        |.........10|
        |0........01|
        |10........0|
        S_{-} is of the form:
        |00000....01|
        |10000.....0|
        |01000.....0|
        |..........0|
        |..........0|
        |.......1000|
        |0.......100|
        |00.......10|'''
    lam = (delt*nu_d)/(delx*delx)
    x_size =  1+int((x_end - x_start)/delx)
    #S_{+}
    S_plus = np.zeros([x_size,x_size])
    S_plus[x_size-1,0] = 1
    for i in range(x_size-1):
        S_plus[i,i+1] = 1
    #S_{-}
    S_minus = np.zeros([x_size,x_size])
    S_minus[0,x_size-1] = 1
    for i in range(1,x_size):
        S_minus[i,i-1] = 1
    I = np.identity(x_size)
    #L matrix
    L = I*(1+lam) - (lam/2)*(S_plus+S_minus)
    #Q matrix
    Q = I*(1-lam) + (lam/2)*(S_plus+S_minus)
    L_inv = np.linalg.inv(L)
    u_n_column = u_n.reshape([x_size,1])
    u_nplus1 = np.matmul(L_inv*Q,u_n_column)
    return u_nplus1.reshape([1,x_size])

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
        u_n = evolution_array[t-1,:]
        u_nplus1 = u_i_nplus1(x_start,x_end,del_t,del_x,nu_d,u_n)
        evolution_array[t,:] = u_nplus1

    u_x_3 = evolution_array[t_size,:]
    
    #error
    glob_err = np.zeros(x_size)
    for x in range(x_size):
        glob_err[x] = E_i_n(u_x_3[x],v_x_3[x]+10e-5)
    #excluding end points coz analytical solution is 0 there
    plt.scatter(x_array[1:-1],glob_err[1:-1])
    plt.title("(i) Global fraction error for del_x = %5.3f, del_t = %5.3f"% (del_x,del_t))
    plt.show()

