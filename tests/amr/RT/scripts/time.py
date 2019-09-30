import matplotlib.pyplot as p
import arepo
import numpy as np

toinch = 0.393700787
p.figure(figsize=np.array([14.7,1.4*14.7])*toinch, dpi=300)

ax1 = p.subplot2grid((17, 5), (0, 0), rowspan=8)
ax2 = p.subplot2grid((17, 5), (0, 1), rowspan=8)
ax3 = p.subplot2grid((17, 5), (0, 2), rowspan=8)
ax4 = p.subplot2grid((17, 5), (0, 3), rowspan=8)
ax5 = p.subplot2grid((17, 5), (0, 4), rowspan=8)

ax6 = p.subplot2grid((17, 5), (8, 0), rowspan=8)
ax7 = p.subplot2grid((17, 5), (8, 1), rowspan=8)
ax8 = p.subplot2grid((17, 5), (8, 2), rowspan=8)
ax9 = p.subplot2grid((17, 5), (8, 3), rowspan=8)
ax10 = p.subplot2grid((17, 5), (8, 4), rowspan=8)

cbar = p.subplot2grid((17,5),(16,0), colspan=4)
cbar2 = p.subplot2grid((17,5),(16,4), colspan=4)

axes = [ax1,ax2,ax3, ax4,ax5,ax6, ax7, ax8, ax9, ax10]
nums = [0,5, 10,15, 20, 25, 30, 35, 40]

cmin = 0.9
cmax = 2.1

for j in np.arange(len(nums)):
    sn = arepo.Simulation("../output/snap_%.3d.hdf5"%nums[j])
    
    sn.plot_AMRslice(sn.rho,axes=axes[j],colorbar=False, cmap = 'parula',res=2048, vmin=cmin, vmax=cmax, gradient=sn.grar)
    axes[j].xaxis.set_visible(False)
    axes[j].yaxis.set_visible(False)
    axes[j].set_frame_on(False)
    axes[j].annotate("t = %.2g"%sn.time, [0.02,1.87],color='black')



arr = np.arange(0,256)
arr = np.array([arr,])
cbar.imshow(arr,aspect='auto',extent=(cmin,cmax,0,1),cmap='parula',rasterized=True, origin="lower")
cbar.set_frame_on(True)
cbar.set_xlabel("$\\rho$")
cbar.set_yticks([])





j = 9
sn = arepo.Simulation("../output/snap_%.3d.hdf5"%40)
    
sn.plot_AMRslice(sn.amrlevel,axes=axes[j],colorbar=False, cmap = 'RdBu',res=2048)
axes[j].xaxis.set_visible(False)
axes[j].yaxis.set_visible(False)
axes[j].set_frame_on(False)
axes[j].annotate("t = %.2g"%sn.time, [0.02,1.87],color='black')


minlevel = sn.amrlevel.min()
maxlevel = sn.amrlevel.max()

arr = np.arange(0,256)
arr[:128] = 0.
arr[128:] = 1.
arr = np.array([arr,])
cbar2.imshow(arr,aspect='auto',extent=(minlevel, maxlevel,0,1),cmap='RdBu',rasterized=True, origin="lower")
cbar2.set_frame_on(True)
cbar2.set_xlabel("AMR level")
cbar2.set_yticks([])
cbar2.set_xticks([minlevel+0.25, maxlevel-0.25])
cbar2.set_xticklabels(["%d"%minlevel, "%d"%maxlevel])

p.tight_layout(h_pad=1., w_pad=1.)
p.savefig("RT_time.pdf", dpi=300)
p.show()
