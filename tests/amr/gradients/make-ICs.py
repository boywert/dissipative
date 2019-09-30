#!/usr/bin/env python

"""
2D test for density gradient computation with AMR
"""

import numpy as np
from collections import namedtuple
from grad import analytic_rho

from gadget import loader

Grid = namedtuple("Grid", "ntot ids x y z vel rho mass vol u")

np.random.seed(42)

gamma = 1.4
boxlength = 1.0
pressure = 1.0

def make_half_grid(N, half):
    ndim = 2
    Ntot = N**ndim / 2
    dx = boxlength / float(N)
    x, y = np.meshgrid(0.5*half + (np.arange(N/2)+0.5)*dx, (np.arange(N)+0.5)*dx, indexing="ij")
    x = x.flatten()
    y = y.flatten()
    z = np.zeros_like(x)
    vel = np.zeros([Ntot,3])
    rho = analytic_rho(x, y)
    vol = (dx**ndim) * np.ones_like(x)
    mass = rho * vol
    ids = np.arange(Ntot, dtype="i")
    u = pressure / ((gamma - 1.0) * rho)
    return Grid(ntot=Ntot, ids=ids, x=x, y=y, z=z, vel=vel, rho=rho, mass=mass, u=u, vol=vol)


level = 6
g1 = make_half_grid(1 << (level+0), 0)
g2 = make_half_grid(1 << (level+1), 1)

def combine(field):
    return np.concatenate([getattr(g1, field), getattr(g2, field)], axis=0)

num_part = np.zeros(6)
num_part[0] = g1.ntot + g2.ntot
num_part[5] = 0

ics = loader.ICs("./ICs.hdf5", num_part)
ics.addField("rho", [1,0,0,0,0,0])
ics.addField("vol", [1,0,0,0,0,0])
ics.part0.pos[:,0] = combine("x")
ics.part0.pos[:,1] = combine("y")
ics.part0.pos[:,2] = combine("z")
ics.part0.vel[:,:] = combine("vel")
ics.part0.rho[:] = combine("rho")
ics.part0.vol[:] = combine("vol")
ics.part0.mass[:] = combine("mass")
ics.part0.id[:] = np.arange(g1.ntot + g2.ntot, dtype="i")
ics.part0.u[:] = combine("u")
ics.write()
