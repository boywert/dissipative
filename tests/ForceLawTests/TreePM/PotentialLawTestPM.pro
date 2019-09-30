
rmin = 0.0001
rmax = 100.0



Basedir = "./output/"
fname = basedir + "forcetest.txt"


spawn,"wc "+fname,result
lines=long(result)
lines=lines(0)


da=dblarr(28, LINES)

openr,1,fname
readf,1,da
close,1


type = da(0,*)
task = da(1,*)
id   = da(2,*)
ti =   da(3,*)


x = da(4,*)
y = da(5,*) 
z = da(6,*) 

r = da(7,*)

; direct summation force
fx = da(8,*)
fy = da(9,*)
fz = da(10,*)

;direct short-range force
fsx = da(11,*)
fsy = da(12,*)
fsz = da(13,*)

;direct long-range force
flx = da(14,*)
fly = da(15,*)
flz = da(16,*)


glx = fx - fsx
gly = fy - fsy
glz = fz - fsz


;tree force
ftx = da(17,*)
fty = da(18,*)
ftz = da(19,*)


;pm force
fpx = da(20,*)
fpy = da(21,*)
fpz = da(22,*)

;direct summation potential
p = da(23,*)
;direct short-range potential
ps = da(24,*)
;direct long-range potential
pl = da(25,*)

;tree potential
pt = da(26,*)

;pm potential
pm = da(27,*)


;treePM total potential
ptm = pt + pm




;TreePM total force
ffx = ftx + fpx
ffy = fty + fpy
ffz = ftz + fpz

f = sqrt(fx^2 + fy^2 + fz^2)
fs = sqrt(fsx^2 + fsy^2 + fsz^2)
fl = sqrt(flx^2 + fly^2 + flz^2)
gl = sqrt(glx^2 + gly^2 + glz^2)

ft = sqrt(ftx^2 + fty^2 + ftz^2)
fp = sqrt(fpx^2 + fpy^2 + fpz^2)
ff = sqrt(ffx^2 + ffy^2 + ffz^2)



window, xsize=1300, ysize=1000
!p.multi=[0,3,2]



plot, r, p, psym=3, /xlog, charsize=2.5, xrange=[rmin, rmax], xstyle=1, title="total potential"
oplot, r, ptm, psym=3, color=255


plot, r, ps, psym=3, /xlog, charsize=2.5, xrange=[rmin, rmax], xstyle=1, title="short-range potential", yrange=[-200,10], ystyle =1 
oplot, r, pt, psym=3, color=255



plot, r, pm, psym=3, /xlog, charsize=2.5, xrange=[rmin, rmax], xstyle=1, title="long-range potential"

oplot, r, pl, psym=3
oplot, r, pm, psym=3, color=255



dp = abs(p-ptm)
plot, r, abs(dp/p), psym=3, /xlog, /ylog, charsize=2, xrange=[rmin, rmax], xstyle=1, title="relative error total potential", yrange=[1.0e-3, max(abs(dp/p))]


dp = abs(ps-pt)
plot, r, abs(dp/ps), psym=3, /xlog, /ylog, charsize=2, xrange=[rmin, rmax], xstyle=1, title="relative error short range potential", yrange=[1.0e-8, max(abs(dp/pt))]



dp = abs(pl-pm)
plot, r, abs(dp/pl), psym=3, /xlog, /ylog, charsize=2, xrange=[rmin, rmax], xstyle=1, title="relative error long range potential", yrange=[1.0e-8, max(abs(dp/pl))]



end
