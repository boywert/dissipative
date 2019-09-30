import numpy as np

import arepo






L_x = 1.
L_y = 1.
N_x = 1024
N_y = 1024

P = 2.5
omega_0 = 0.1
sigma = 0.05/np.sqrt(2)
GAMMA = 5./3.

dx = L_x/N_x
dy = L_y/N_y


sn = arepo.ICs("ics_static.hdf5", [N_x*N_y,0,0,0,0,0], precision=np.float64)

mesh = np.meshgrid(np.arange(N_x),np.arange(N_y))

posx = (mesh[0]*dx).reshape(N_x*N_y)+0.5*dx
posy = (mesh[1]*dy).reshape(N_x*N_y)+0.5*dy

mesh = np.meshgrid(np.arange(N_x),np.arange(N_y))
posx2 = (mesh[0]*dx).reshape(N_x*N_y)+0.5*dx
posy2 = (mesh[1]*dy).reshape(N_x*N_y)+0.5*dy

sn.part0.pos[:N_x*N_y/4,0] = posx[:N_x*N_y/4]
sn.part0.pos[:N_x*N_y/4,1] = posy[:N_x*N_y/4]
sn.part0.mass[:N_x*N_y/4] = 1.*dx*dy
sn.part0.u[:N_x*N_y/4] = P/((GAMMA-1)*1.)
sn.part0.vel[:N_x*N_y/4,0] = -0.5

sn.part0.pos[N_x*N_y/4:N_x*N_y/4+N_x*N_y/2,0] = posx2[N_x*N_y/4:N_x*N_y*3/4]
sn.part0.pos[N_x*N_y/4:N_x*N_y/4+N_x*N_y/2,1] = posy2[N_x*N_y/4:N_x*N_y*3/4]
sn.part0.mass[N_x*N_y/4:N_x*N_y/4+N_x*N_y/2] = 2.*dx*dy
sn.part0.u[N_x*N_y/4:N_x*N_y/4+N_x*N_y/2] = P/((GAMMA-1)*2.)
sn.part0.vel[N_x*N_y/4:N_x*N_y/4+N_x*N_y/2,0] = +0.5

sn.part0.pos[-N_x*N_y/4:,0] = posx[-N_x*N_y/4:]
sn.part0.pos[-N_x*N_y/4:,1] = posy[-N_x*N_y/4:]
sn.part0.mass[-N_x*N_y/4:] = 1.*dx*dy
sn.part0.u[-N_x*N_y/4:] = P/((GAMMA-1)*1.)
sn.part0.vel[-N_x*N_y/4:,0] = -0.5

sn.part0.pos[:,2] = 0.


sn.part0.vel[:,1] = omega_0 * np.sin(4*np.pi*sn.part0.pos[:,0])*(np.exp(-(sn.part0.pos[:,1]-0.25)**2*0.5/(sigma**2))+np.exp(-(sn.part0.pos[:,1]-0.75)**2*0.5/(sigma**2)))
sn.part0.vel[:,2] = 0.

sn.part0.id[:] = np.arange(N_x*N_y)+1



sn.write()

