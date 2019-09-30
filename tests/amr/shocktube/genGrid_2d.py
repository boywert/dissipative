import numpy as np

import arepo



L_x = 2./4
L_y = 2.
N_x = 1024/4
N_y = 1024


GAMMA = 1.4

dx = L_x/N_x
dy = L_y/N_y


sn = arepo.ICs("ics.hdf5", [N_x*N_y,0,0,0,0,0], precision=np.float64)

mesh = np.meshgrid(np.arange(N_x),np.arange(N_y))

posx = (mesh[0]*dx).reshape(N_x*N_y)+0.5*dx
posy = (mesh[1]*dy).reshape(N_x*N_y)+0.5*dy

sn.part0.pos[:,0] = posx
sn.part0.pos[:,1] = posy

rho = np.zeros((N_x*N_y))
rho[:N_y/2*N_x] = 1.
rho[N_y/2*N_x:] = 0.125

sn.part0.mass[:] = rho*dx*dy

sn.part0.u[:N_y/2*N_x] = 1. /((GAMMA-1)* 1. )
sn.part0.u[N_y/2*N_x:] = 0.1 /((GAMMA-1)* 0.125 )


sn.part0.vel[:N_y/2*N_x,1] = 0.75
sn.part0.vel[N_y/2*N_x:,1] = 0.

sn.part0.id[:] = np.arange(N_x*N_y)+1


sn.header.boxsize = L_y

sn.write()

