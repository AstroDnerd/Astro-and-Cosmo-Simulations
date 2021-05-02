""" @package examples/wave_2d/create.py
Code that creates 2d wave test problem;

created by Nikhil Bisht, 23.3.21; Last modified: 8.4.21
"""

""" load libraries """
import sys    # system specific calls
import numpy as np    # scientific computing package
import h5py    # hdf5 format

simulation_directory = str(sys.argv[1])
print("examples/wave_2d/create.py: creating ICs in directory " + simulation_directory)

""" initial condition parameters """
FilePath = simulation_directory + '/IC.hdf5'

## ensure calculations happen with predefined precision
FloatType = np.float64  # double precision: np.float64, for single use np.float32
IntType = np.int32 # integer type

## computational domain
Boxsize = FloatType(1.0)
if len(sys.argv) > 3:
    CellsPerDimension = IntType(sys.argv[3])
else:
    CellsPerDimension = IntType(40)
NumberOfCells = CellsPerDimension * CellsPerDimension

## initial state
density_0 = FloatType(1.0)
velocity_radial_0 = FloatType(0.2)	#Radial outflow
velocity_theta_0 = FloatType(0.0)
pressure_0 = FloatType(3.0) / FloatType(5.0)
gamma = FloatType(5.0) / FloatType(3.0)
gamma_minus_one = gamma - FloatType(1.0)
delta = FloatType(1e-6)    # relative density perturbation
uthermal_0 = pressure_0 / density_0 / gamma_minus_one

""" set up grid: cartesian 2d grid """
## spacing
dx = Boxsize / FloatType(CellsPerDimension)
## position of first and last cell
pos_first, pos_last = 0.5 * dx, Boxsize - 0.5 * dx

## set up grid
Grid1d = np.linspace(pos_first, pos_last, CellsPerDimension, dtype=FloatType)
xx, yy = np.meshgrid(Grid1d, Grid1d)
Pos = np.zeros([NumberOfCells, 3], dtype=FloatType)
Pos[:,0] = xx.reshape(NumberOfCells)
Pos[:,1] = yy.reshape(NumberOfCells)
## calculate distance from center
xPosFromCenter = (Pos[:,0] - 0.5 * Boxsize)
yPosFromCenter = (Pos[:,1] - 0.5 * Boxsize)
Radius = np.sqrt( xPosFromCenter**2 + yPosFromCenter**2 )
Volume = np.full(NumberOfCells, dx, dtype=FloatType)

""" set up hydrodynamical quantitites """
## set up unperturbed system; density, velocity and specific internal energy
Density = np.full(Pos.shape[0], density_0, dtype=FloatType)
Velocity = np.zeros(Pos.shape, dtype=FloatType)
Uthermal = np.full(Pos.shape[0], uthermal_0, dtype=FloatType)

## perturbations
Density *=  FloatType(1.0) + delta * np.sin( FloatType(2.0) * FloatType(np.pi) * Pos[:,0] / Boxsize )
Velocity[:,0] = velocity_radial_0 * xPosFromCenter / Radius
Velocity[:,1] = velocity_radial_0 * yPosFromCenter / Radius
Uthermal *= (Density / density_0)**gamma_minus_one

## mass instead of density needed for input
Mass = Density * Volume

""" write *.hdf5 file; minimum number of fields required by Arepo """
IC = h5py.File(FilePath, 'w')    # open/create file

## create hdf5 groups
header = IC.create_group("Header")    # create header group
part0 = IC.create_group("PartType0")    # create particle group for gas cells

## write header entries
NumPart = np.array([NumberOfCells, 0, 0, 0, 0, 0], dtype=IntType)
header.attrs.create("NumPart_ThisFile", NumPart)
header.attrs.create("NumPart_Total", NumPart)
header.attrs.create("NumPart_Total_HighWord", np.zeros(6, dtype=IntType) )
header.attrs.create("MassTable", np.zeros(6, dtype=IntType) )
header.attrs.create("Time", 0.0)
header.attrs.create("Redshift", 0.0)
header.attrs.create("BoxSize", Boxsize)
header.attrs.create("NumFilesPerSnapshot", 1)
header.attrs.create("Omega0", 0.0)
header.attrs.create("OmegaB", 0.0)
header.attrs.create("OmegaLambda", 0.0)
header.attrs.create("HubbleParam", 1.0)
header.attrs.create("Flag_Sfr", 0)
header.attrs.create("Flag_Cooling", 0)
header.attrs.create("Flag_StellarAge", 0)
header.attrs.create("Flag_Metals", 0)
header.attrs.create("Flag_Feedback", 0)
if Pos.dtype == np.float64:
    header.attrs.create("Flag_DoublePrecision", 1)
else:
    header.attrs.create("Flag_DoublePrecision", 0)

## write cell data
part0.create_dataset("ParticleIDs", data=np.arange(1, NumberOfCells+1) )
part0.create_dataset("Coordinates", data=Pos)
part0.create_dataset("Masses", data=Mass)
part0.create_dataset("Velocities", data=Velocity)
part0.create_dataset("InternalEnergy", data=Uthermal)

## close file
IC.close()

""" normal exit """
sys.exit(0)
