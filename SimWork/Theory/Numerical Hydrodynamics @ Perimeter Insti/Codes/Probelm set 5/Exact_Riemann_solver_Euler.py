#Exact Riemann solver for 1D time dependent Euler Equations for ideal gas

import numpy as np 
import matplotlib.pyplot as plt
import time

#Initial conditions and constants
gamma = 1.4
rho_l = 1.0
u_l = 0.0
p_l = 1.0
rho_r = 0.125
u_r = 0.0
p_r = 0.1
c_l = np.sqrt(gamma*p_l/rho_l)
c_r = np.sqrt(gamma*p_r/rho_r)

'''we get p_star by finding the roots of the following function:
f(p,w_l,w_r) := f_l(p,w_l) + f_r(p,w_r) + u_r - u_l
f_l(p,w_l) = (p-p_l)*sqrt(A_l/(p+B_l)) for p>p_l
           = (2c_l/(gamma-1))*((p/p_l)^((gamma-1)/2gamma)-1) for p<=p_l

and f_r is the same as f_l with r replacing l
and
A_l = 2/((gamma+1)rho_l); B_l = ((gamma-1)/(gamma+1))*p_l
and A_r and B_r are same as A_l and B_l replacing r with l

We use Newton Raphson to find the root with an initial
guess of 0.5(p_l+p_r) and a tolerance of 10^-6'''

A_l = 2/((gamma+1)*rho_l)
B_l = ((gamma-1)/(gamma+1))*p_l
A_r = 2/((gamma+1)*rho_r)
B_r = ((gamma-1)/(gamma+1))*p_r

def f_lr(p,p_rl,rho_rl,c_rl,A_rl,B_rl,gamma):
	if p>p_rl:
		return (p-p_rl)*np.sqrt(A_rl/(p+B_rl))
	else:
		return (2*c_rl/(gamma-1))*((p/p_rl)**((gamma-1)/2*gamma)-1)

def fprime_lr(p,p_rl,rho_rl,c_rl,A_rl,B_rl,gamma):
	if p>p_rl:
		return ((p+2*B_rl+p_rl)/(2*(p+B_rl)))*np.sqrt(A_rl/(p+B_rl))
	else:
		return (c_rl/(gamma*p_rl))*((p_rl/p)**((gamma+1)/2*gamma))
def f_main(p,p_l,rho_l,c_l,A_l,B_l,p_r,rho_r,c_r,A_r,B_r,gamma,u_r,u_l):
	f_l = f_lr(p,p_l,rho_l,c_l,A_l,B_l,gamma)
	f_r = f_lr(p,p_r,rho_r,c_r,A_r,B_r,gamma)
	return f_l+f_r+u_r-u_l
def fprime_main(p,p_l,rho_l,c_l,A_l,B_l,p_r,rho_r,c_r,A_r,B_r,gamma,u_r,u_l):
	fp_l = fprime_lr(p,p_l,rho_l,c_l,A_l,B_l,gamma)
	fp_r = fprime_lr(p,p_r,rho_r,c_r,A_r,B_r,gamma)
	return fp_l+fp_r

def Newton_Raphson_step(p_guess,p_l,rho_l,c_l,A_l,B_l,p_r,rho_r,c_r,A_r,B_r,gamma,u_r,u_l):
	f_main_p = f_main(p_guess,p_l,rho_l,c_l,A_l,B_l,p_r,rho_r,c_r,A_r,B_r,gamma,u_r,u_l)
	f_prime_main_p = fprime_main(p_guess,p_l,rho_l,c_l,A_l,B_l,p_r,rho_r,c_r,A_r,B_r,gamma,u_r,u_l)
	return (p_guess - (f_main_p/f_prime_main_p))

def get_p_star(p_l,rho_l,c_l,A_l,B_l,p_r,rho_r,c_r,A_r,B_r,gamma,u_r,u_l):
	p_guess = 0.5*(p_l+p_r)
	p_new = Newton_Raphson_step(p_guess,p_l,rho_l,c_l,A_l,B_l,p_r,rho_r,c_r,A_r,B_r,gamma,u_r,u_l)
	while np.abs(p_new-p_guess)>1e-06:
		p_guess = p_new
		p_new = Newton_Raphson_step(p_guess,p_l,rho_l,c_l,A_l,B_l,p_r,rho_r,c_r,A_r,B_r,gamma,u_r,u_l)

	return p_new
'''We can get u_star from p_star using:
u_star = 0.5*(u_l+u_r) + 0.5*(f_r(p_star)-f_p(p_star))'''
def get_u_star(p_star,p_l,rho_l,c_l,A_l,B_l,p_r,rho_r,c_r,A_r,B_r,gamma,u_r,u_l):
	f_l = f_lr(p_star,p_l,rho_l,c_l,A_l,B_l,gamma)
	f_r = f_lr(p_star,p_r,rho_r,c_r,A_r,B_r,gamma)
	return (0.5*(u_l+u_r) + 0.5*(f_r-f_l))

p_star_val = get_p_star(p_l,rho_l,c_l,A_l,B_l,p_r,rho_r,c_r,A_r,B_r,gamma,u_r,u_l)
u_star_val = get_u_star(p_star_val,p_l,rho_l,c_l,A_l,B_l,p_r,rho_r,c_r,A_r,B_r,gamma,u_r,u_l)
#some extra quantities

#left case
S_l = u_l-(c_l*np.sqrt(((gamma+1)/(2*gamma))*(p_star_val/p_l) + (gamma-1)/(2*gamma)))
rho_star_l_case11 = rho_l*((p_star_val/p_l)+(gamma-1)/(gamma+1))/((p_star_val/p_l)*((gamma-1)/(gamma+1))+1)
rho_star_l_case12 = rho_l*((p_star_val/p_l)**(1/gamma))
c_star_l = c_l*((p_star_val/p_l)**((gamma-1)/(2*gamma)))
S_Hl = u_l-c_l
S_Tl = u_star_val-c_star_l

#right case
S_r = u_r+(c_r*np.sqrt(((gamma+1)/(2*gamma))*(p_star_val/p_r) + (gamma-1)/(2*gamma)))
rho_star_r_case11 = rho_r*((p_star_val/p_r)+(gamma-1)/(gamma+1))/((p_star_val/p_r)*((gamma-1)/(gamma+1))+1)
rho_star_r_case12 = rho_r*((p_star_val/p_r)**(1/gamma))
c_star_r = c_r*((p_star_val/p_r)**((gamma-1)/(2*gamma)))
S_Hr = u_r+c_r
S_Tr = u_star_val+c_star_r

#solving for the primitives
def get_primitive(S,u_star,p_star):
	if S<=u_star:
		#left side
		if p_star>p_l:
			#left shock
			if (S>=S_l):
				return (rho_star_l_case11,u_star,p_star)
			else:
				return (rho_l,u_l,p_l)
		else:
			#left rarefaction wave
			if S<=S_Hl:
				return (rho_l,u_l,p_l)
			elif (S>S_Hl and S<=S_Tl):
				rho_fan = rho_l*(((2/(gamma+1))+((gamma-1)/((gamma+1)*c_l))*(u_l-S))**(2/(gamma-1)))
				u_fan = (2/(gamma+1))*(c_l+(gamma-1)/2 + S)
				p_fan = p_l*(((2/(gamma+1))+((gamma-1)/((gamma+1)*c_l))*(u_l-S))**(2*gamma/(gamma-1)))
				return (rho_fan,u_fan,p_fan)
			else:
				return (rho_star_l_case12,u_star,p_star)
	else:
		#right side
		if p_star>p_r:
			#right shock
			if (S<=S_r):
				return (rho_star_r_case11,u_star,p_star)
			else:
				return (rho_r,u_r,p_r)
		else:
			#right rarefaction wave
			if S<=S_Tr:
				return (rho_star_r_case12,u_star,p_star)
			elif (S>S_Tr and S<=S_Hr):
				rho_fan = rho_r*(((2/(gamma+1))-((gamma-1)/((gamma+1)*c_r))*(u_r-S))**(2/(gamma-1)))
				u_fan = (2/(gamma+1))*(-c_r+(gamma-1)/2 + S)
				p_fan = p_r*(((2/(gamma+1))-((gamma-1)/((gamma+1)*c_r))*(u_r-S))**(2*gamma/(gamma-1)))
				return (rho_fan,u_fan,p_fan)
			else:
				return (rho_r,u_r,p_r)

#Solving the Euler Equation
t_initial = 0.01
delta_t = 0.01
t_max = 0.5
gridsize = 101
t_size = int((t_max-t_initial)/delta_t)+1
t_grid = np.zeros(t_size)
x_grid = np.linspace(-1,1,num=gridsize)
rho = np.zeros([t_size,gridsize])
v = np.zeros([t_size,gridsize])
p = np.zeros([t_size,gridsize])
e = np.zeros([t_size,gridsize])
t_cnt = 0
t=t_initial
while t<t_max:
	x_cnt=0
	for x in x_grid:
		(rho[t_cnt,x_cnt],v[t_cnt,x_cnt],p[t_cnt,x_cnt]) = get_primitive(x/t,u_star_val,p_star_val)
		e[t_cnt,x_cnt] = (p[t_cnt,x_cnt]/(gamma-1))+0.5*rho[t_cnt,x_cnt]*v[t_cnt,x_cnt]**2
		x_cnt+=1
	t_grid[t_cnt] = t
	t+=delta_t
	t_cnt+=1

#plots
fig, axs = plt.subplots(4,figsize=(8,8))
for t in t_grid:
	t_index = int((t-t_initial)/delta_t)
	
	fig.suptitle('Plot of Primitives at time t=%0.2f'%(t))
	axs[0].plot(x_grid, rho[t_index,:],'-o')
	axs[0].set_ylabel("Density")
	axs[1].plot(x_grid, v[t_index,:],'-o')
	axs[1].set_ylabel("Velocity")
	axs[2].plot(x_grid, p[t_index,:],'-o')
	axs[2].set_ylabel("Pressure")
	axs[3].plot(x_grid, e[t_index,:],'-o')
	axs[3].set_ylabel("Energy")
	axs[3].set_xlabel("x Grid")
	plt.draw()
	plt.pause(0.01)
	time.sleep(0.1)
	axs[0].clear()
	axs[1].clear()
	axs[2].clear()
	axs[3].clear()
	#plt.clf()
