PicName = "frame"
PicOutName = "pic"

BaseDir ="../picdat/"

PicDir  ="../pics/"

flag = 0

;for Num = 0,4095 do begin
for Num = 0,678,678 do begin


  exts='0000'
  exts=exts+strcompress(string(num),/remove_all)
  exts=strmid(exts,strlen(exts)-4,4)

  fname=BaseDir +PicName+"_" +exts +".dat"
  fname=strcompress(fname,/remove_all)
  
  openr,1,fname

  PixelsX=0L
  PixelsY=0L
  readu,1, PixelsX, PixelsY

min_massmap=0.0
max_massmap=0.0
min_veldispmap=0.0
max_veldispmap=0.0
min_densmap=0.0
max_densmap=0.0
min_velzmap=0.0 
max_velzmap=0.0


  readu,1, min_massmap, max_massmap
  readu,1, min_veldispmap, max_veldispmap
  readu,1, min_densmap, max_densmap
  readu,1, min_velzmap, max_velzmap

;goto,ende
  
  MassMap= fltarr(PixelsY, PixelsX)
  VelDispMap= fltarr(PixelsY, PixelsX)
  DensMap= fltarr(PixelsY, PixelsX)
  VelZMap= fltarr(PixelsY, PixelsX)
  readu,1,MassMap
  readu,1,VelDispMap
  readu,1,DensMap
  readu,1,VelZMap
  close,1



print, "Min:", min_massmap, min(MassMap)
print, "Max:", max_massmap, max(MassMap)



  
  MassMap= transpose(MassMap)
  VelDispMap= transpose(VelDispMap)
  DensMap= transpose(DensMap)
  VelZMap= transpose(VelZMap)
  
  
  Value = DensMap
  
  ma1=Max(Value) ;/20
  mi1=Min(Value) ;*5
  
  print, "mi1=", mi1
  print, "ma1=", ma1
  
  
  ind=where(Value gt ma1)
  if ind(0) ne -1 then begin
     Value(ind)=ma1
  endif
  ind=where(Value lt mi1)
  if ind(0) ne -1 then begin
     Value(ind)=mi1
  endif
  
  imageDens= byte((alog(Value)-alog(mi1))/(alog(ma1)-alog(mi1))* 255.9)
  
  Value = VelDispMap / DensMap

  ma2= max(Value)
  mi2= min(Value)
  
 ; mi2 = 5.0
 ; ma2 = 500.0
  
  print, "mi = ", mi2, " km/sec"
  print, "ma = ", ma2, " km/sec"
  
  
  ind=where(Value gt ma2)
  if ind(0) ne -1 then begin
     Value(ind)=ma2
  endif
  ind=where(Value lt mi2)
  if ind(0) ne -1 then begin
     Value(ind)=mi2
  endif
  
  imageVelDisp= long((alog(Value)-alog(mi2))/(alog(ma2)-alog(mi2))* 359.9)
  

;goto,ende

  cols=256
  ColPlane=bytarr(cols,360,3)
  
  openr,1,"colcube_6-10_saturated.dat"
  readu,1,ColPlane
  close,1
  

  JpegPic=bytarr(PixelsX,PixelsY,3)
  
  for i=0, PixelsX-1 do begin
     for j=0, PixelsY-1 do begin
        JPegPic(i,j,0)= ColPlane(imageDens(i,j),imageVelDisp(i,j),0)
        JPegPic(i,j,1)= ColPlane(imageDens(i,j),imageVelDisp(i,j),1)
        JPegPic(i,j,2)= ColPlane(imageDens(i,j),imageVelDisp(i,j),2)
     endfor
  endfor


  exts='0000'
  exts=exts+strcompress(string(num),/remove_all)
  exts=strmid(exts,strlen(exts)-4,4)
  
  
  fname = picdir + "/"+picoutname+"_"+exts+".jpg"
  
  write_jpeg, fname, JPegPic, true=3, quality=95
  
  print, fname
  
  if flag eq 0 then begin
 ;    window, xsize=PixelsX,ysize=pixelsY
     flag = 1
  endif
;  tv, JPegPic, true=3

  
endfor

ende:
end
