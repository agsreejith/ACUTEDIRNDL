; NAME:
;      COSMIC
;
; PURPOSE:
;      Returns the 2d image array after putting in n number of rndom cosmic ray hits 
; CALLING SEQUENCE:
;      Result = COSMIC(number,count)
;
; INPUTS:
;      number = Number of cosmic rays in a frame.
;      count  = Count of cosmic ray.
;      nx,ny  = number of piuxels in x and y
;
; OUTPUT:
;      Returns 2d image array in counts/s/cm2
; REQUIRES:
;      The number of cosmic ray hits to be added.
;
; PROCEDURE:
;      Creates a 2d image array with random cosmic ray hits
; MODIFICATION HISTORY:
;      created 04.12.2017 by A. G. Sreejith
;

function cosmic,number,count,nx,ny
;the length of cosmic ray is fixed as 6 pixel

width=6
im=dblarr(nx,ny)

x_val=(nx-1)*randomu(seed,number)
y_val=(ny-1)*randomu(seed,number)
direction=360*randomu(seed,number)
direction=0.01745329252d*direction ;convertion to radians
for i=0,number-1 do begin
  x=fix(x_val[i])
  y=fix(y_val[i])
  im(x,y)=count
  for j=0, 5 do begin
  x_shft=round(j*cos(direction[i]))
  y_shft=round(j*sin(direction[i]))  
  x_n=x+x_shft
  y_n=y+y_shft
  ;boundary condition to limit consmic ray within the image
  if x_n gt nx-1 then x_n = nx-1
  if y_n gt ny-1 then y_n = ny-1
  if x_n lt 0 then x_n = 0
  if y_n lt 0 then y_n = 0
  im[x_n,y_n]=count
  endfor
endfor
return,im
end



