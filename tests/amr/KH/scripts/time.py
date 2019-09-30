import matplotlib.pyplot as p
import arepo
import numpy as np

toinch = 0.393700787
p.figure(figsize=np.array([14.7,5])*toinch, dpi=300)

ax1 = p.subplot2grid((1, 3), (0, 0))
ax2 = p.subplot2grid((1, 3), (0, 1))
ax3 = p.subplot2grid((1, 3), (0, 2))


axes = [ax1,ax2,ax3]
nums = [0,100,200]

for j in np.arange(len(nums)):
    sn = arepo.Simulation("../output_static/snap_%.3d.hdf5"%nums[j])
    
    sn.plot_AMRslice(sn.rho,axes=axes[j],colorbar=False)
    axes[j].xaxis.set_visible(False)
    axes[j].yaxis.set_visible(False)
    axes[j].set_frame_on(False)
    axes[j].annotate("t = %.2g"%sn.time, [0.02,0.9],color='white')

p.tight_layout()
p.show()
p.savefig("KH_time.pdf", dpi=300)
