; NAME:
;      COSMIC
;
; PURPOSE:
;      Returns the 2d image array after putting in n number of rndom cosmic ray hits 
; CALLING SEQUENCE:
;      Result = COSMIC(number)
;
; INPUTS:
;      number = Number of cosmic rays in a frame.
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

function cosmic, number
;the length of cosmic ray is fixed as 6 pixel

width=6
im=dblarr(2048,515)

x_val=2048*randomu(seed,number)
y_val=515*randomu(seed,number)
direction=360*randomu(seed,number)
direction=0.01745329252d*direction ;convertion to radians
for i=0,number-1 do begin
  x=fix(x_val[i])
  y=fix(y_val[i])
  im(x,y)=1000000
  for j=0, 5 do begin
  x_shft=round(j*cos(direction[i]))
  y_shft=round(j*sin(direction[i]))  
  x_n=x+x_shft
  y_n=y+y_shft
  im(x_n,y_n)=1000000
  endfor
endfor
return,im
end



