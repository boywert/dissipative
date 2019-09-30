import numpy as np

import arepo



L_x = 0.5
L_y = 2.
N_x = 256
N_y = 4*N_x

P = 2.5
g = -0.1
rho1 = 1.
rho2 = 2.

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
rho[:N_y/2*N_x] = rho1
rho[N_y/2*N_x:] = rho2

sn.part0.mass[:] = rho*dx*dy

sn.part0.u[:] = (P + g * (sn.pos[:,1] - 1.) * rho)/((GAMMA-1)* rho )
sn.part0.id[:] = np.arange(N_x*N_y)+1


sn.part0.vel[:,1] = 0.0025 * (1. - np.cos(4.*np.pi*sn.part0.pos[:,0])) * (1. - np.cos(np.pi*sn.part0.pos[:,1]))

sn.header.boxsize = 1.

sn.write()

