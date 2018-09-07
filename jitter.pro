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

function jitter, image_array, dx, dy, dz, exp_time,pscale 

  if not keyword_defined(image_array ) then image_array = dblarr(2048,515)
  if not keyword_defined(dx          ) then dx          = 0
  if not keyword_defined(dy          ) then dy          = 0
  if not keyword_defined(exp_time    ) then exp_time    = 1

s=size(image_array)
nx=s[1]
ny=s[2]

rand=fix(20*randomu(seed,1))
ran_val=rand[0]
for rd=0,ran_val do count_x=randomn(seed,exp_time)
rand=fix(20*randomu(seed,1))
ran_val=rand[0]
for rd=0,ran_val do count_y=randomn(seed,exp_time)
rand=fix(20*randomu(seed,1))
ran_val=rand[0]
for rd=0,ran_val do angle=randomn(seed,exp_time)
dx=dx/pscale ;converting to per pixel
dy=dy/pscale ;converting to per pixel
angle=angle*dz/3600.
x=fix(count_x*dx/2.35)
y=fix(count_y*dy/2.35)
image_out=dblarr(nx,ny)
for i=0, exp_time-1 do begin
im_new=shift(image_array,x[i],y[i])  
im_rot = ROT(im_new, angle[i], 1, /INTERP, /MISSING)
image_out=image_out+im_rot
endfor

return,image_out
end