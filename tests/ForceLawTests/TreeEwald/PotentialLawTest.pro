
rmin = 0.0001
rmax = 100.0



Basedir = "./output/"
fname = basedir + "forcetest.txt"


spawn,"wc "+fname,result
lines=long(result)
lines=lines(0)


da=dblarr(16, LINES)

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

;tree force
ftx = da(11,*)
fty = da(12,*)
ftz = da(13,*)



;direct summation potential
p = da(14,*)
;direct tree potential
pt = da(15,*)



f = sqrt(fx^2 + fy^2 + fz^2)
ft = sqrt(ftx^2 + fty^2 + ftz^2)



window, xsize=1300, ysize=800
!p.multi=[0,2,1]



plot, r, p, psym=3, /xlog, charsize=2.5, xrange=[rmin, rmax], xstyle=1, title="total potential"
oplot, r, pt, psym=3, color=255


dp = abs(p-pt)
plot, r, abs(dp/p), psym=3, /xlog, /ylog, charsize=2, xrange=[rmin, rmax], xstyle=1, title="relative error total potential", yrange=[1.0e-3, max(abs(dp/p))]


end
