; NAME:
;      CUTE_SPECSPREAD
;
; PURPOSE:
;      Returns the 2d image array after spreading the 1d spectrum
; CALLING SEQUENCE:
;      Result = CUTE_SPECSPREAD(file_pos, file_spread, counts, nx, ny)
;
; INPUTS:
;      file_pos = The file specifying positional information of spectrum asuming a linear function inputs from the file are y intersept, slope and start pixel on CCD.
;      file_spread = The file specifying the spread in y axis (cross dispersion) the spread is currently assumed as guassian.
;      ccd_count = The counts/pixel of the spectrum, same number of element as pixels in dispersion direction.
;      nx,ny = The number of pixels in x and y. 
;      s_pos = Position on the slit.
;      
; OUTPUT:
;      Returns 2d image arrayin counts/s/cm2
; REQUIRES:
;      The position input and spread in ccd.
;
; PROCEDURE:
;      Convert 1d spectra into 2d spectra on the detector
; MODIFICATION HISTORY:
;      created 01.12.2017 by A. G. Sreejith
;      

function keyword_defined, key

  if n_elements(key) eq 0 then return, 0 else return, 1

end


function cute_specspread, file_pos, file_spread, ccd_count, nx, ny,s_pos

if not keyword_defined(file_pos   ) then file_pos    = 'extra\position.txt'
if not keyword_defined(file_spread) then file_spread = 'extra\cuteback.txt'
if not keyword_defined(nx         ) then nx          = 2048
if not keyword_defined(ny         ) then ny          = 515

length=file_lines(file_pos)
data_pos=dblarr(3,length)
openr,1,file_pos
readf,1,data_pos
close,1
  ;default central position
  intercept=-0.561
  slope=-4.34e-4
  start= 257

if s_pos eq 0 then begin
    intercept=data_pos[0,1]
    slope=data_pos[0,2]
    start=data_pos[0,0]
endif

if s_pos eq 4 then begin
    intercept=data_pos[1,1]
    slope=data_pos[1,2]
    start=data_pos[0,0]
endif

if s_pos eq 8 then begin
    intercept=data_pos[2,1]
    slope=data_pos[2,2]
    start=data_pos[0,0]
endif
length=file_lines(file_spread)
openr,1,file_spread
data_spread=dblarr(5,length)
readf,1,data_spread
close,1
x_spread=data_spread[3,*]
xrms=2*data_spread[2,*]
;transalation in pixels to center of the detector default 257
im=dblarr(nx,ny)

for x=0,nx-1 do begin
  y=slope*x+intercept
  y=fix(y)
  im[x,y+start]=ccd_count(x)
  npix=fix(x_spread[x])
  width=fix(xrms[x])
  xspread=psf_Gaussian( NPIXEL=npix, FWHM=width , /DOUBLE, /NORMALIZE);spread in cross dispesion direction guassian
  for j=0,npix-1 do begin
    im[x,(y+start+j-(npix/2))]=xspread[j]*ccd_count(x)
  endfor
endfor

return,im
end
