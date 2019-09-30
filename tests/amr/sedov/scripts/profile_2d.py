import matplotlib.pyplot as p
import arepo
import numpy as np

toinch = 0.393700787
p.figure(figsize=np.array([14.7,10])*toinch, dpi=300)



sn = arepo.Simulation("../output_2d/snap_%.3d.hdf5"%10)

sn.plot_radprof(sn.rho, label="t = %.2g"%sn.time, bins = 200)

t = np.loadtxt("output_2d.dat",skiprows=2)

p.plot(t[:,1],t[:,2], label="analytical solution",zorder=0)

p.legend(loc=4, frameon=False, fontsize=9)
p.tight_layout()
p.show()
p.savefig("sedov_profile_2d.pdf", dpi=300)
