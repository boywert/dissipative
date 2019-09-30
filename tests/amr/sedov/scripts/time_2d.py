import matplotlib.pyplot as p
import arepo
import numpy as np
from mpl_toolkits.axes_grid1 import AxesGrid

toinch = 0.393700787
"""
p.figure(figsize=np.array([14,14])*toinch, dpi=300)

ax1 = p.subplot2grid((3, 3), (0, 0))
ax2 = p.subplot2grid((3, 3), (0, 1))
ax3 = p.subplot2grid((3, 3), (0, 2))

ax4 = p.subplot2grid((3, 3), (1, 0))
ax5 = p.subplot2grid((3, 3), (1, 1))
ax6 = p.subplot2grid((3, 3), (1, 2))

ax7 = p.subplot2grid((3, 3), (2, 0))
ax8 = p.subplot2grid((3, 3), (2, 1))
ax9 = p.subplot2grid((3, 3), (2, 2))
"""

fig = p.figure(figsize =np.array([14.7,14.7*1.3])*toinch)
axes = AxesGrid(fig, 111,
                    nrows_ncols = (3, 3),
                    axes_pad = 0.1,
                    share_all=False,
                    #label_mode = "L",
                    cbar_location = "bottom",
                    cbar_mode="single",
                    )

#axes = [ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8,ax9]
nums = np.array([3,9,12,15,18,21,24,27,30])

cmin = 0.
cmax = 5.

for j in np.arange(len(nums)):
    sn = arepo.Simulation("../output_2d/snap_%.3d.hdf5"%nums[j])
    
    sn.plot_AMRslice(sn.rho,axes=axes[j],colorbar=False,vmin=cmin,vmax=cmax, cmap='parula')
    axes[j].xaxis.set_visible(False)
    axes[j].yaxis.set_visible(False)
    axes[j].set_frame_on(False)
    axes[j].annotate("t = %.2g"%sn.time, [0.02,0.9],color='white')


arr = np.arange(0,256)
arr = np.array([arr,])
axes.cbar_axes[0].imshow(arr,aspect='auto',extent=(cmin,cmax,0,1),cmap='parula',rasterized=True, origin="lower")
axes.cbar_axes[0].set_frame_on(True)
axes.cbar_axes[0].set_xlabel("$\\rho$")


p.tight_layout()
p.show()
p.savefig("sedov_time_2d.pdf", dpi=300)
