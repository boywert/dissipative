import matplotlib.pyplot as p
import arepo
import numpy as np

toinch = 0.393700787
p.figure(figsize=np.array([12,16])*toinch, dpi=300)

ax1 = p.subplot2grid((5, 1), (0, 0))
ax2 = p.subplot2grid((5, 1), (1, 0))
ax3 = p.subplot2grid((5, 1), (2, 0))
ax4 = p.subplot2grid((5, 1), (3, 0))


exact = np.loadtxt("e1rpex.out")




axes = [ax1,ax2,ax3,ax4]
sims = ['../output/snap_010.hdf5','../output_arepo/snap_010.hdf5','../output_albada/snap_010.hdf5','../output_minbee/snap_010.hdf5','../output_superbee/snap_010.hdf5','../output_vanleer/snap_010.hdf5']
label = ["None", "Arepo", "Albada", "Minbee", "Superbee", "Vanleer"]

sims = ['../output_superbee/snap_020.hdf5' ,'../output_vanleer/snap_020.hdf5','../output_albada/snap_020.hdf5','../output_minbee/snap_020.hdf5','../output_arepo/snap_020.hdf5']
label = [ "\\small{SUPERBEE}", "van Leer","van Albada", "\\small{MINBEE}", "Springel 2010"]




lines = []

l, = ax1.plot(exact[:,0],exact[:,1],ls='--', lw=0.5)
lines.append(l)
label.insert(0, 'analytical solution')

ax2.plot(exact[:,0],exact[:,3],ls='--', lw=0.5)
ax3.plot(exact[:,0],exact[:,2],ls='--', lw=0.5)
ax4.plot(exact[:,0],exact[:,4],ls='--', lw=0.5)

for j in np.arange(len(sims)):
    sn = arepo.Simulation(sims[j])
    sn.set_center([0.25,1.2])
    g = sn.get_AMRline(sn.rho, gradient=sn.grar,axis=1, box=1.)
    l, = ax1.plot(g["x2"]-0.7,g['grid'], lw=0.5)
    lines.append(l)
    g = sn.get_AMRline(sn.pres, gradient=sn.grap, axis=1, box=1.)
    ax2.plot(g["x2"]-0.7,g['grid'],lw = 0.5)
    g = sn.get_AMRline(sn.vel_y, gradient = sn.grav[:,3:6], axis=1, box=1.)
    ax3.plot(g["x2"]-0.7,g['grid'],lw = 0.5)
    g = sn.get_AMRline(sn.u, axis=1, box=1.)
    ax4.plot(g["x2"]-0.7,g['grid'],lw = 0.5)
    
ax1.set_ylim(0.,1.1)
ax2.set_ylim(0.,1.1)
ax3.set_ylim(-0.1,1.5)
ax4.set_ylim(1.9,3.7)


ax1.set_xlim(0.45,0.6)
ax2.set_xlim(0.45,0.6)
ax3.set_xlim(0.45,0.6)
ax4.set_xlim(0.45,0.6)

    
p.figlegend(lines,label,"lower center", ncol = 2)
p.tight_layout()
p.show()
p.savefig("shocktube_compare.pdf", dpi=300)
