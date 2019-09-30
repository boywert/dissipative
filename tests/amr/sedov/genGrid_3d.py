import numpy as np

import arepo



L_x = 1.
L_y = 1.
L_z = 1.
N_x = 64
N_y = 64
N_z = 64

P = 1e-4

GAMMA = 5./3.

dx = L_x/N_x
dy = L_y/N_y
dz = L_z/N_z

sn = arepo.ICs("ics_3d.hdf5", [N_x*N_y*N_z,0,0,0,0,0], precision=np.float64)

mesh = np.mgrid[0:N_x, 0:N_y, 0:N_z]

posx = (mesh[0]*dx).reshape(N_x*N_y*N_z)+0.5*dx
posy = (mesh[1]*dy).reshape(N_x*N_y*N_z)+0.5*dy
posz = (mesh[2]*dz).reshape(N_x*N_y*N_z)+0.5*dz

sn.part0.pos[:,0] = posx
sn.part0.pos[:,1] = posy
sn.part0.pos[:,2] = posz
sn.part0.mass[:] = 1.*dx*dy*dz
sn.part0.u[:] = P/((GAMMA-1)*1.)


E0 = 1./sn.part0.mass[0]/8.

sn.part0.u[(N_x/2-1)*N_y*N_z + (N_y/2-1)*N_z + (N_z/2-1)] = E0
sn.part0.u[(N_x/2-1)*N_y*N_z + (N_y/2-1)*N_z + (N_z/2)] = E0
sn.part0.u[(N_x/2-1)*N_y*N_z + (N_y/2)*N_z + (N_z/2-1)] = E0
sn.part0.u[(N_x/2-1)*N_y*N_z + (N_y/2)*N_z + (N_z/2)] = E0
sn.part0.u[(N_x/2)*N_y*N_z + (N_y/2-1)*N_z + (N_z/2-1)] = E0
sn.part0.u[(N_x/2)*N_y*N_z + (N_y/2-1)*N_z + (N_z/2)] = E0
sn.part0.u[(N_x/2)*N_y*N_z + (N_y/2)*N_z + (N_z/2-1)] = E0
sn.part0.u[(N_x/2)*N_y*N_z + (N_y/2)*N_z + (N_z/2)] = E0



sn.part0.id[:] = np.arange(N_x*N_y*N_z)+1

sn.header.boxsize = 1.

sn.write()

