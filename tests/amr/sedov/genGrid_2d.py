import numpy as np

import arepo



L_x = 1.
L_y = 1.
N_x = 128
N_y = 128

P = 1e-4

GAMMA = 5./3.

dx = L_x/N_x
dy = L_y/N_y


sn = arepo.ICs("ics_2d.hdf5", [N_x*N_y,0,0,0,0,0], precision=np.float64)

mesh = np.meshgrid(np.arange(N_x),np.arange(N_y))

posx = (mesh[0]*dx).reshape(N_x*N_y)+0.5*dx
posy = (mesh[1]*dy).reshape(N_x*N_y)+0.5*dy

sn.part0.pos[:,0] = posx
sn.part0.pos[:,1] = posy
sn.part0.mass[:] = 1.*dx*dy
sn.part0.u[:] = P/((GAMMA-1)*1.)


E0 = 1./sn.part0.mass[0]/4.

sn.part0.u[(N_x/2-1)*N_y + (N_y/2-1)] = E0
sn.part0.u[(N_x/2-1)*N_y + (N_y/2)] = E0
sn.part0.u[(N_x/2)*N_y + (N_y/2-1)] = E0
sn.part0.u[(N_x/2)*N_y + (N_y/2)] = E0


sn.part0.pos[:,2] = 0.
sn.part0.vel[:,1] = 0.
sn.part0.vel[:,0] = 0.


sn.part0.id[:] = np.arange(N_x*N_y)+1

sn.header.boxsize = 1.

sn.write()

