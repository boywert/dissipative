
Base= "../../../output/"

DETAILED_TIMING_GRAVWALK    =     0
DETAILED_TIMING_STELLARDENSITY =  1


Ncpu = 240
Step = 131315
Part = DETAILED_TIMING_STELLARDENSITY



window, xsize=1200, ysize=1000

for t=0, Ncpu-1 do begin
   
   fname = Base + "/timings_detailed_" + strcompress(string(t),/remove_all) + ".txt"


   spawn,"wc "+fname,result
   lines=long(result)
   lines=lines(0)
   
   
   en=fltarr(6,LINES)
   
   
   openr,1,fname
   readf,1,en
   close,1
   
   tistep = en(0,*)
   tibin = en(1,*)
   codepart = en(2,*)
   mode = en(3,*)

   tstart = en(4,*)
   tend = en(5,*)

   if t eq 0 then begin

      ind = where((tistep eq STEP) and (codepart eq PART) and (mode eq 2))
      
      timebin = max(tibin(ind))

      ind = where((tistep eq STEP) and (tibin eq timebin) and (codepart eq PART) and (mode eq 2))

      tt = tend(ind(0))

      plot, [0],[0], psym=3, xrange=[0,tt*1.05], yrange =[0,NCPU], xstyle=1, ystyle=1

      oplot, tt*[1,1], [0,NCPU], color=255*256L^2, thick=2.0

   endif

   ind = where((tistep eq STEP) and (tibin eq timebin) and (codepart eq PART) and (mode eq 0))
   if ind(0) ne -1 then begin
      for rep =0,n_elements(ind)-1 do begin
         polyfill, [tstart(ind(rep)),tend(ind(rep)), tend(ind(rep)), tstart(ind(rep)), tstart(ind(rep))], [t,t, t+0.3, t+0.3, t], color=255
      endfor
   endif

   ind = where((tistep eq STEP) and (tibin eq timebin) and (codepart eq PART) and (mode eq 1))
   if ind(0) ne -1 then begin
      for rep =0,n_elements(ind)-1 do begin
         polyfill, [tstart(ind(rep)),tend(ind(rep)), tend(ind(rep)), tstart(ind(rep)), tstart(ind(rep))], [t,t, t+0.3, t+0.3, t], color=255*256L^1+255
      endfor
   endif

endfor



end

