import matplotlib.pyplot as p
import arepo
import numpy as np

toinch = 0.393700787
p.figure(figsize=np.array([14.7*2./3,5])*toinch, dpi=300)

ax1 = p.subplot2grid((1, 2), (0, 0))
ax2 = p.subplot2grid((1, 2), (0, 1))



axes = [ax1,ax2]
sims = ['../output_amr/snap_200.hdf5','../output_amr_nosmoothing/snap_200.hdf5']
label = ["AMR", "AMR no smoothing"]

for j in np.arange(len(sims)):
    sn = arepo.Simulation(sims[j])
    
    sn.plot_AMRslice(sn.amrlevel,axes=axes[j],colorbar=False, cmap="RdBu")
    axes[j].xaxis.set_visible(False)
    axes[j].yaxis.set_visible(False)
    axes[j].set_frame_on(False)
    axes[j].annotate(label[j], [0.02,0.9],color='white')

p.tight_layout()
p.show()
p.savefig("KH_amrlevel.pdf", dpi=300)
