# import libraries
import numpy as np
from netCDF4 import Dataset
from scipy import interpolate

# load from .dat file
forwardmodel = np.loadtxt('Data/bestfit_Amatrice.dat')
forwardmodel = forwardmodel[2:] # first two are dummies

# using 13 (along strike) x 8 (along dip) sparse grid
T0I = forwardmodel[0:104]
fsI = forwardmodel[104:208]
DcI = forwardmodel[208:312]

T0I_grid = np.zeros((13,8))
fsI_grid = np.zeros((13,8))
DcI_grid = np.zeros((13,8))
for i in range(13):
    T0I_grid[i] = (np.flipud((T0I[i::13]).T)).T
    fsI_grid[i] = (np.flipud((fsI[i::13]).T)).T
    DcI_grid[i] = (np.flipud((DcI[i::13]).T)).T
    
# bilinear interpolate T0I, DcI and fsI
x = np.linspace(0, 30000, 13)
z = np.linspace(0, 14000, 8)
kind = 'linear'
f_T0I = interpolate.interp2d(x, z, T0I_grid.T, kind=kind)
f_DcI = interpolate.interp2d(x, z, DcI_grid.T, kind=kind)
f_fsI = interpolate.interp2d(x, z, fsI_grid.T, kind=kind)
grid = 100
xx = np.arange(0, 30000+grid, grid)
zz = np.arange(0, 14000+grid, grid) # using fault top depth at 0 m
T0I_grid_new = (f_T0I(xx, zz)).T # pre-stress (Tau_i)
DcI_grid_new = (f_DcI(xx, zz)).T # critical distance (Dc)
fsI_grid_new = (f_fsI(xx, zz)).T # static - dynamic friction difference

# set parameters
dip = 45 
sI = 8.520 * np.sin(dip*np.pi/180) # normal stress depth gradient 8.520 MPa/km
normalstress1d = np.linspace(zz[0]*sI/1000*1e6,zz[-1]*sI/1000*1e6,zz.shape[0])
normalstress1d[normalstress1d<=0.1e6] = 0.1e6

T0_n = T0I_grid_new * 0
for k in range (0,T0_n.shape[0]):
    T0_n[k] = normalstress1d

fdI = 0.4 # mu_d in SeisSol (dynamic fric. coeff.)
cohesion = -0.5e6 # cohesion in SeisSol
T0_l = T0I_grid_new * 0
T0_m = T0I_grid_new + T0_n * fdI # pre-stress + dynamic friction stress

# dynamic rupture parameters
fs = fsI_grid_new + fdI # static friction coeff. (fsI = static - dynamic fric. coeff)
mud = fs * 0.0 + fdI # dynamic friction coeff.
d_c = DcI_grid_new * 1 # critical distance
coh = fs * cohesion # cohesion

# set 2d fault
nx = 301
nz = 141

tn0 = -T0_n.T
ts0 = T0_l.T
td0 = T0_m.T
mus0 = fs.T
mud0 = mud.T
d_c0 = d_c.T
coh0 = cohesion.T

xmin = -18000
xmax = 12000
zmin = 0
zmax = -14000

x = np.linspace(xmin,xmax,nx)
z = np.linspace(zmin,zmax,nz)

# create nc file T_n, T_s, T_d
fout = Dataset('Amatrice_fault_DR_2d_TnTsTd.nc','w',format='NETCDF4')
fout.createDimension('x',nx)
fout.createDimension('z',nz)
fout.createVariable('x','f4',('x',))
fout.createVariable('z','f4',('z',))
fout.variables['x'][:] = x;
fout.variables['z'][:] = z;

fault = np.dtype([('T_n','f4'),('T_s','f4'),('T_d','f4')])
data = fout.createCompoundType(fault,'fault')
vv = fout.createVariable('data',data,dimensions=('z','x'))

ii=0
for k in range(0,nz):
    for i in range(0,nx):
        vv[k,i] = (tn0[k,i],ts0[k,i],td0[k,i])
        ii += 1
fout.close()

# create nc file mu_s, mu_d, d_c, cohesion
fout = Dataset('Amatrice_fault_DR_2d_MusMudDcCoh.nc','w',format='NETCDF4')
fout.createDimension('x',nx)
fout.createDimension('z',nz)
fout.createVariable('x','f4',('x',))
fout.createVariable('z','f4',('z',))
fout.variables['x'][:] = x;
fout.variables['z'][:] = z;

fault = np.dtype([('mu_s','f4'),('mu_d','f4'),('d_c','f4'),('cohesion','f4')])
data = fout.createCompoundType(fault,'fault')
vv = fout.createVariable('data',data,dimensions=('z','x'))

ii=0
for k in range(0,nz):
    for i in range(0,nx):
        vv[k,i] = (mus0[k,i],mud0[k,i],d_c0[k,i],coh0[k,i])
        ii += 1
fout.close()
