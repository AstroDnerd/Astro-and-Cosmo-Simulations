""" @package ./examples/wave_2d/check.py
Code that checks results of 2d wave propagation problem

created by Nikhil Bisht, 23.03.2021 (custom example)
"""

""" load libraries """
import sys    # system specific calls
import numpy as np    # scientific computing package
import h5py    # hdf5 format
import os     # file specific calls
import matplotlib.pyplot as plt    # needs to be active for plotting!
plt.rcParams['text.usetex'] = True

simulation_directory = str(sys.argv[1])
print("wave_2d: checking simulation output in directory " + simulation_directory) 

FloatType = np.float64  # double precision: np.float64, for single use np.float32
IntType = np.int32 # integer type

makeplots = True
if len(sys.argv) > 2:
  if sys.argv[2] == "True":
    makeplots = True
  else:
    makeplots = False

""" open initial conditiions to get parameters """
try:
    data = h5py.File(simulation_directory + "/IC.hdf5", "r")
except:
    print("could not open initial  conditions!")
    exit(-1)
Boxsize = FloatType(data["Header"].attrs["BoxSize"])
NumberOfCells = IntType(data["Header"].attrs["NumPart_Total"][0]) 
CellsPerDimension = np.sqrt(NumberOfCells) ## 2d sim

""" maximum L1 error after one propagation; empirically motivated """
DeltaMaxAllowed = 5e-5 * FloatType(NumberOfCells)**-2

""" initial state -- copied from create.py """
density_0 = FloatType(1.0)
velocity_radial_0 = FloatType(0.2)	#Radial outflow
velocity_theta_0 = FloatType(0.0)
pressure_0 = FloatType(3.0) / FloatType(5.0)
gamma = FloatType(5.0) / FloatType(3.0)
gamma_minus_one = gamma - FloatType(1.0)
delta = FloatType(1e-6)    # relative density perturbation
uthermal_0 = pressure_0 / density_0 / gamma_minus_one


""" 
    loop over all output files; need to be at times when analytic
    solution equals the initial conditions
"""
i_file = 0
status = 0
error_data = []
while True:
    """ try to read in snapshot """
    directory = simulation_directory+"/output/"
    filename = "snap_%03d.hdf5" % (i_file)
    try:
        data = h5py.File(directory+filename, "r")
    except:
        break
    
    """ get simulation data """
    ## simulation data
    time = FloatType(data["Header"].attrs["Time"])
    ## simulation data
    Pos = np.array(data["PartType0"]["CenterOfMass"], dtype = FloatType)
    VoronoiPos = np.array(data["PartType0"]["Coordinates"], dtype = FloatType)
    Density = np.array(data["PartType0"]["Density"], dtype = FloatType)
    Mass = np.array(data["PartType0"]["Masses"], dtype = FloatType)
    Velocity = np.array(data["PartType0"]["Velocities"], dtype = FloatType)
    Uthermal = np.array(data["PartType0"]["InternalEnergy"], dtype = FloatType)
    CellVolume = Mass / Density
    
    xPosFromCenter = (Pos[:,0] - 0.5 * Boxsize)
    yPosFromCenter = (Pos[:,1] - 0.5 * Boxsize)
    Radius = np.sqrt( xPosFromCenter**2 + yPosFromCenter**2 )
    
    vRad = Velocity[:,0] * xPosFromCenter / Radius + Velocity[:,1] * yPosFromCenter / Radius
    
    """ calculate analytic solution at new cell positions """
    Density_ref = np.full(Pos.shape[0], density_0, dtype=FloatType)
    Velocity_ref = np.zeros(Pos.shape, dtype=FloatType)
    Uthermal_ref = np.full(Pos.shape[0], uthermal_0, dtype=FloatType)
    ## perturbations
    Density_ref *=  FloatType(1.0) + delta * np.sin( FloatType(2.0) * FloatType(np.pi) * Pos[:,0] / Boxsize )
    Velocity_ref[:,0] = velocity_radial_0
    Uthermal_ref *= (Density / density_0)**gamma_minus_one

    """ compare data """
    ## density
    abs_delta_dens = np.abs(Density - Density_ref) / Density_ref
    L1_dens = np.average(abs_delta_dens, weights=CellVolume)
    
    ## velocity, here, use absolute error (velocity should be zero! check only x-vel, taking all components dilutes the norm!)
    abs_delta_vel = np.abs(Velocity - Velocity_ref)[:,0]
    L1_vel = np.average(abs_delta_vel, weights=CellVolume)
    
    ## internal energy
    abs_delta_utherm = np.abs(Uthermal - Uthermal_ref) / Uthermal_ref
    L1_utherm = np.average(abs_delta_utherm, weights=CellVolume)

    """ printing results """
    print("wave_1d: L1 error of " + filename +":")
    print("\t density: %g" % L1_dens)
    print("\t velocity: %g" % L1_vel)
    print("\t specific internal energy: %g" % L1_utherm)
    print("\t tolerance: %g for %d cells" % (DeltaMaxAllowed, NumberOfCells) )
    
    error_data.append(np.array([L1_dens, L1_vel, L1_utherm], dtype=FloatType))
    
    
    """ criteria for failing the test """
    if L1_dens > DeltaMaxAllowed or L1_vel > DeltaMaxAllowed or L1_utherm > DeltaMaxAllowed:
        status = 1
    
    if makeplots and i_file > 0:
      if not os.path.exists( simulation_directory+"/plots" ):
        os.mkdir( simulation_directory+"/plots" )
      ## plot:
      fig = plt.figure( figsize=np.array([7.0,3.5]), dpi=300 )
      
      Nplot = 256
      from scipy import spatial # needed for KDTree that we use for nearest neighbour search
      Edges1d = np.linspace(0., Boxsize, Nplot+1, endpoint=True, dtype=FloatType)
      Grid1d = 0.5 * (Edges1d[1:] + Edges1d[:-1])
      xx, yy = np.meshgrid(Grid1d, Grid1d)
      Grid2D = np.array( [xx.reshape(Nplot**2), yy.reshape(Nplot**2)] ).T
      
      vor = spatial.Voronoi( VoronoiPos[:,:2] )
      dist, cells = spatial.KDTree( VoronoiPos[:,:2] ).query( Grid2D, k=1 )
      
      ax  = plt.axes( [0.05,0.15,0.35,0.7] )
      pc  = ax.pcolormesh( Edges1d, Edges1d, Velocity[cells,0].reshape((Nplot,Nplot)), rasterized=True )
      cax = plt.axes( [0.42,0.15,0.01,0.7] )
      plt.colorbar( pc, cax=cax )
      ax.set_title( 'Radial velocity' )
      spatial.voronoi_plot_2d( vor, ax=ax, line_colors='w', show_points=False, show_vertices=False, line_width=0.5 )
      ax.set_xlim( 0, Boxsize )
      ax.set_ylim( 0, Boxsize )
      
      error = Velocity[cells,0]-Velocity_ref[cells,0]
      ax  = plt.axes( [0.53,0.15,0.35,0.7] )
      pc  = ax.pcolormesh( Edges1d, Edges1d, error.reshape((Nplot,Nplot)), rasterized=True )
      cax = plt.axes( [0.90,0.15,0.01,0.7] )
      plt.colorbar( pc, cax=cax )
      ax.set_title( 'Radial velocity error' )
      spatial.voronoi_plot_2d( vor, ax=ax, line_colors='w', show_points=False, show_vertices=False, line_width=0.5 )
      ax.set_xlim( 0, Boxsize )
      ax.set_ylim( 0, Boxsize )
      
      plt.text( 0.5, 0.92, "$t=%4.1f,\ L_1=%5.1e$" % (time,L1_dens), ha='center', size=12, transform=fig.transFigure )

      if not os.path.exists( simulation_directory+"/plots" ):
        os.mkdir( simulation_directory+"/plots" )

      print(simulation_directory+"/plots/figure_%03d.pdf" % (i_file) )
      fig.savefig(simulation_directory+"/plots/figure_%03d.pdf" % (i_file), dpi=300)
      plt.close(fig)
    i_file += 1
np.savetxt(simulation_directory+"/error_%d.txt"%NumberOfCells, np.array(error_data))

""" normal exit """
sys.exit(status) 
