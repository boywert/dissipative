
time_dat = dblarr(2, 4096)
openr, 1, "ExpansionList_4096"
readf,1,time_dat
close,1





ClusX=  29.6306      
ClusY= 212.360      
ClusZ= 447.576


BoxSize = 500.0

Depth =    50.0

Angle1  =  60.0 ; full opening angle of field of view 
Angle2 =   20.0

PixelsX =    1920
PixelsY =    1080



R=  0.05
xc= 0.5
yc= 0.5

L = 1.0 + 2*!PI*R
N = 2500



T = indgen(N)/float(N-1)*L

X=fltarr(N)
Y=fltarr(N)

ind = where(T lt xc)


d = (xc - T(ind))

X(ind) = T(ind)
Y(ind) = tanh((2*T(ind)/xc-1)*4.0)*R/2 +r/2 + yc

;plot, X(ind), Y(ind)


ind = where((T ge xc) and (T lt (xc + 2*!PI*r)))
phi = (T(ind) -xc)/r

X(ind) = xc + r* sin(phi)
Y(ind) = yc + r *cos(phi)




ind = where((T ge (xc + 2*!PI*r)))

X(ind) = T(ind) - (xc + 2*!PI*r)  + xc
Y(ind) = tanh((2*(T(ind)-(xc + 2*!PI*r))/xc-1)*4.0)*(-R/2) +r/2 + yc



;plot, X, Y


KX=fltarr(N)
KY=fltarr(N)
Kx(*)=1

;ind = where((T lt (xc + 0.75*2*!PI*r)))
ind = where((T lt 10000.0))

vx = xc-X(ind)
vy = yc-Y(ind)

vv=sqrt(vx^2+vy^2)
Kx(ind) = vx/vv
Ky(ind) = vy/vv


kkx =kx
kky =ky




kkx  =  tanh((2*(T-1.05)/0.1)-1)  

kky  =  sqrt(1-kkx^2)



kx += kkx+1
ky += kky


;vx = X(ind(1:*))-X(ind(1:*)-1)
;vy = Y(ind(1:*))-y(ind(1:*)-1)
;vv=sqrt(vx^2+vy^2)
;Kx(ind(1:*)) = vx/vv
;Ky(ind(1:*)) = vy/vv


; plot, t, ky
;oplot, t, kx, color=255




slowfac = 0.8

sigma = 2.5*!pi*r


K = sqrt(!PI/2) * sigma * (errorf( (L-L/2)/(Sqrt(2)*sigma)) + errorf( (L/2)/(Sqrt(2)*sigma)))



AA = L /(L-slowfac*K)
BB = slowfac*AA

TT = AA*T - BB*sqrt(!PI/2) * sigma * (errorf( (T-L/2)/(Sqrt(2)*sigma)) + errorf( (L/2)/(Sqrt(2)*sigma)))


;plot, t, tt
print, max(t), max(tt)






t0 = 0.5 + R*!PI/2.0
t1 = 0.5 + R*3*!PI/2.0

wt =  R*!PI/8.0


AngleTab  =  angle1 + (angle2-angle1)*0.5 * (tanh((T-t0)/wt) + 1) + (angle1-angle2)*0.5 * (tanh((T-t1)/wt) + 1)
 



common cosmicdata, Omega, OmegaLambda


Omega=0.25D
OmegaLambda = 0.75D
Hubble=0.73

a=1.0 
tage=h0t0(a)

print,"a= ",a
print,"z= ",1/a-1
print,"t= ",tage

Tage=tage/0.1 *0.98/Hubble
print, "T= ", Tage, "Gyr"

NN= 10000
alist = (indgen(NN) + 1.0d) / NN
tlist = dblarr(NN)
for i=0L, NN-1 do begin
  tlist(i) = h0t0(alist(i))
endfor

tbase = interpol(tlist, alist, 1/(50.0+1))
tmax = max(tlist)

PICS = 4096

ascalepic=dblarr(PICS)

openw, 1, "piclist.txt", width = 1000


for i=0, PICS-1 do begin

    ti = (L/PICS)*(i+0.5)

    ascale = interpol(alist, tlist, tbase + ((tmax-tbase) * i)/(PICS-1))

    time_dat(1,i)  = (i+0.5)/PICS
    ascalepic(i)= ascale

    ti = interpol(TT, T, ti)


    Angle = interpol(AngleTab, T, ti)

    cameraX = interpol(X, T, ti) 
    cameraY = interpol(Y, T, ti) 
    cameraZ = 0.5

    cameraX *= BoxSize
    cameraY *= BoxSize
    cameraZ *= BoxSize


    cameraX += (ClusX - 250.0)
    cameraY += (ClusY - 250.0)
    cameraZ += (ClusZ - 250.0)



 ;   LengthX =  2 * Depth * tan(0.5* Angle/180.0 *!PI)
 ;   LengthY =  LengthX * float(PixelsY)/PixelsX

     LengthY =  0.75 * 2 * Depth * tan(0.5* Angle/180.0 *!PI)
     LengthX =  LengthY * float(PixelsX)/PixelsY


    dirX = interpol(KX, T, ti) 
    dirY = interpol(KY, T, ti) 
    dirZ = 0.0

    upX = 0.0
    upY = 0.0
    upZ = 1.0

    rightX = dirY*upZ - dirZ*upY
    rightY = dirZ*upX - dirX*upZ
    rightZ = dirX*upY - dirY*upX
        
 

    printf,1, i, ascalepic(i),  $ 
       cameraX, cameraY, cameraZ, $
           -rightX, -rightY, -rightZ, $
           upX, upY, upZ, $
           Depth, Angle

endfor

 
close,1


;plot, time_dat(1,*), time_dat(0,*)
;oplot, time_dat(1,*), ascalepic, color=255





ende:
end
