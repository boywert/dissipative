
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

;pm pot
ppm = da(23,*)
pew = da(24,*)

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



plot, r, f, psym=3, /xlog, /ylog, charsize=2.5, xrange=[rmin, rmax], xstyle=1, title="total force"
oplot, r, ff, psym=3, color=255

plot, r, fs, psym=3, /xlog, /ylog, charsize=2.5, xrange=[rmin, rmax], xstyle=1, title="short-range force", yrange=[max(fs)/1.0e6, max(fs)]
oplot, r, ft, psym=3, color=255



plot, r, fl, psym=3, /xlog, /ylog, charsize=2.5, xrange=[rmin, rmax], xstyle=1, title="long-range force"
oplot, r, fp, psym=3, color=255

oplot, r, gl, psym=3, color=255*256L+255



df = sqrt((fx-ffx)^2 + (fy-ffy)^2 + (fz-ffz)^2)
plot, r, df/f, psym=3, /xlog, /ylog, charsize=2, xrange=[rmin, rmax], xstyle=1, title="relative error total force", yrange=[1.0e-3, max(df/f)]


df = sqrt((fsx-ftx)^2 + (fsy-fty)^2 + (fsz-ftz)^2)
plot, r, df/fs, psym=3, /xlog, /ylog, charsize=2, xrange=[rmin, rmax], xstyle=1, title="relative error short range force", yrange=[1.0e-8, max(df/fs)]



df = sqrt((flx-fpx)^2 + (fly-fpy)^2 + (flz-fpz)^2)
plot, r, df/fp, psym=3, /xlog, /ylog, charsize=2, xrange=[rmin, rmax], xstyle=1, title="relative error long rangeforce", yrange=[1.0e-8, max(df/fp)]



end
