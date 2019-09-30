import matplotlib.pyplot as p
import arepo
import numpy as np
from mpl_toolkits.axes_grid1 import AxesGrid

toinch = 0.393700787
"""
p.figure(figsize=np.array([14.7,5])*toinch, dpi=300)

ax1 = p.subplot2grid((1, 3), (0, 0))
ax2 = p.subplot2grid((1, 3), (0, 1))
ax3 = p.subplot2grid((1, 3), (0, 2))
"""

fig = p.figure(figsize =np.array([14.7,10.])*toinch)
axes = AxesGrid(fig, 111,
                    nrows_ncols = (1, 2),
                    axes_pad = 0.1,
                    share_all=False,
                    #label_mode = "L",
                    cbar_location = "bottom",
                    cbar_mode="single",
                    )




sims = ['../output_static/snap_200.hdf5', '../output_static_grad/snap_200.hdf5']
label = ["Superbee", "Springel 2010"]

cmin = 0.9
cmax = 2.1

for j in np.arange(len(sims)):
    sn = arepo.Simulation(sims[j])
    
    sn.plot_AMRslice(sn.rho,axes=axes[j],res=2048, gradient=sn.grar, cmap="parula", vmin=cmin, vmax=cmax, colorbar=False)
    axes[j].xaxis.set_visible(False)
    axes[j].yaxis.set_visible(False)
    axes[j].set_frame_on(False)
    axes[j].annotate(label[j], [0.02,0.9],color='white')

arr = np.arange(0,256)
arr = np.array([arr,])
axes.cbar_axes[0].imshow(arr,aspect='auto',extent=(cmin,cmax,0,1),cmap='parula',rasterized=True, origin="lower")
axes.cbar_axes[0].set_frame_on(True)
axes.cbar_axes[0].set_xlabel("$\\rho$")
    
p.tight_layout()

p.savefig("KH_compare_highres.pdf", dpi=300)
p.show()
