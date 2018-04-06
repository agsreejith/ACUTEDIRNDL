; NAME:
;      JITTER
;
; PURPOSE:
;      Returns the 2d image array after simulating the effect of jitter of the spacecraft specified
; CALLING SEQUENCE:
;      Result = JITTER(in_image, dx, dy, exptime)
;
; INPUTS:
;      in_image = Input image array to which jitter effect has to be added.
;      dx,dy,dz = The RMS shift of pixels in x and y and in z (angular) in arcseconds.
;      exptime = Exposure time of the image.
; OUTPUT:
;      Returns 2d image array in counts/s/cm2
; REQUIRES:
;      The RMS shift in X and Y and the exposure time. 
;
; PROCEDURE:
;      Silulates the effect of jitter.
; MODIFICATION HISTORY:
;      created 01.12.2017 by A. G. Sreejith
;      modified 27.12.2017 by A. G. Sreejith
;       
;

function jitter, image_array, dx, dy, dz, exp_time 

  if not keyword_defined(image_array ) then image_array = dblarr(2048,515)
  if not keyword_defined(dx          ) then dx          = 0
  if not keyword_defined(dy          ) then dy          = 0
  if not keyword_defined(exp_time    ) then exp_time    = 1

count_x=randomn(seed,exp_time)
count_y=randomn(seed,exp_time)
angle=randomn(seed,exp_time)
dx=dx/2.5 ;converting to per pixel
dy=dy/2.5 ;converting to per pixel
angle=angle*dz/3600
x=fix(count_x*dx/2.35)
y=fix(count_y*dy/2.35)
image_out=dblarr(2048,515)
for i=0, exp_time-1 do begin
im_new=shift(image_array,x[i],y[i])  
im_rot = ROT(im_new, angle[i], 1, /INTERP, /MISSING)
image_out=image_out+im_rot
endfor
return,image_out
end