import numpy as np
import math

class Integrator:
    def __init__(self, xMin, xMax, N):
        self.xmax = xMax
        self.xmin = xMin
        self.n = N
        self.integ=0
    
            
    def integrate(self):
        x = np.linspace(self.xmin,self.xmax,self.n)
        delx = (self.xmax-self.xmin)/self.n
        fx = (x**2)*np.exp(-x)*np.sin(x)
        self.integ = np.sum(fx)*delx
        
        
    def show(self):
        print(self.integ)

        

examp = Integrator(1,3,200000)
examp.integrate()
examp.show()
