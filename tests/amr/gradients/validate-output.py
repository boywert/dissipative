#!/usr/bin/env python

import numpy as np
import sys
from argparse import ArgumentParser
import matplotlib.pyplot as plt

import arepo
from grad import analytic_grad_rho

parser = ArgumentParser()
parser.add_argument("snapshot")
args = parser.parse_args()

out = arepo.Simulation(args.snapshot)

# Discard border cells
x, y, z = out.pos.T
margin = 0.1
inside = (x > margin) & (x < 1-margin) & (y > margin) & (y < 1-margin)

theta = np.deg2rad(57)
gx, gy, gz = out.grar[inside,:].T

# Reference gradient (analytic)
grefx, grefy, grefz = analytic_grad_rho(x[inside], y[inside])

err = np.sqrt( (gx-grefx)**2 + (gy-grefy)**2 + (gz-grefz)**2 )

print "Gradient error (max, L2):", np.max(err), np.sqrt(np.mean(err**2))

plt.figure(figsize=(12, 12), dpi=120)
plt.axvline(0.5)
plt.quiver(x[inside], y[inside], gx-np.cos(theta), gy-np.sin(theta), scale=5)
plt.title(args.snapshot)
plt.show()

