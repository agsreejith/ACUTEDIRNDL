; NAME:
;      BADPIX
;
; PURPOSE:
;      Returns the 2d image array after putting in bad and hot pixels including bad/hot coloums/rows
; CALLING SEQUENCE:
;      Result = COSMIC(file)
;
; INPUTS:
;      file = Parameter file for bad/hot pixel positions.
;
; OUTPUT:
;      Returns 2d image array in counts/s/cm2
; REQUIRES:
;      The parameter file for bad/hot pixel.
;
; PROCEDURE:
;      Creates a 2d image array with bad/hot pixels.
; MODIFICATION HISTORY:
;      created 04.12.2017 by A. G. Sreejith
;
function badpix,file
length=file_lines(file)
data=intarr(4,length)
openr,1,file
readf,1,data
close,1
im=dblarr(2048,515)
print,data
for i=0,length-1 do begin
 x=data(0,i)
 y=data(1,i)
 hb=data(2,i)  
 rc=data(3,i)
 print,'rc',rc
 case rc of
   0: if (hb eq 0) then im(x,y) = -1 else if (hb eq 1) then im(x,y) = 10000 else print,'Wrong input in file: check badpixel map file)  
   1: if (hb eq 0) then im(x,*) = -1 else if (hb eq 1) then im(x,*) = 10000 else print,'Wrong input in file: check badpixel map file)
   2: if (hb eq 0) then im(*,y) = -1 else if (hb eq 1) then im(*,y) = 10000 else print,'Wrong input in file: check badpixel map file)
 endcase
endfor
return,im
end
