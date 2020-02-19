; NAME:
;      CUTE_SPECSPREAD
;
; PURPOSE:
;      Returns the 2d image array after spreading the 1d spectrum, Options for spread based on  
; CALLING SEQUENCE:
;      Result = CUTE_SPECSPREAD(file_pos, file_spread, counts, nx, ny)
;
; INPUTS:
;      file_pos = The file specifying positional information of spectrum asuming a linear function inputs from the file are y intersept, slope and start pixel on CCD.
;      file_spread = The file specifying the spread in y axis (cross dispersion) the spread is currently assumed as guassian.
;      ccd_count = The counts/pixel of the spectrum, same number of element as pixels in dispersion direction.
;      nx,ny = The number of pixels in x and y. 
;      s_pos = Position on the slit in arc minutes from center.
;      x_shift = Shift in x in arc minutes from central pixel (mainly for background stars)
;      y_shift = Shift in y in arc minutes from central pixel (mainly for background stars)
; OUTPUT:
;      Returns 2d image arrayin counts/s/cm2
; REQUIRES:
;      The position input and spread in ccd.
;
; PROCEDURE:
;      Convert 1d spectra into 2d spectra on the detector
; MODIFICATION HISTORY:
;      created 01.12.2017 by A. G. Sreejith
;      modified 27.12.2017 by A. G. Sreejith
;       added xdispersion based on actual values from Zemax model. 
;      

function cute_specspread, file_pos, file_spread, ccd_count, nx, ny,s_pos, x_shift,y_shift,scale
close,/all
if not keyword_defined(file_pos    ) then file_pos    = 'extra\position.txt'
if not keyword_defined(file_spread ) then file_spread = 'extra\cuteback.txt'
if not keyword_defined(nx          ) then nx          = 2048
if not keyword_defined(ny          ) then ny          = 515
if not keyword_defined(x_shift     ) then x_shift     = 0
if not keyword_defined(y_shift     ) then y_shift     = 0
if not keyword_defined(scale       ) then scale       = 2.5

n_pix=60.0/double(scale)
length=file_lines(file_pos)
data_pos=dblarr(3,length)
openr,1,file_pos
readf,1,data_pos
close,1
  ;default central position
  intercept=-0.561
  slope=-4.34e-4
  start= 257

if (s_pos ge -2 && s_pos le 2) then begin
    intercept=data_pos[0,1]
    slope=data_pos[0,2]
    start=data_pos[0,0]+fix(s_pos*n_pix)
endif

if ((s_pos gt 2 && s_pos le 8) or (s_pos lt -2 && s_pos ge -8)) then begin
    intercept=data_pos[1,1]
    slope=data_pos[1,2]
    start=data_pos[0,0]+fix(s_pos*n_pix)
endif

if ((s_pos gt 8) or (s_pos lt -8)) then begin
    intercept=data_pos[2,1]
    slope=data_pos[2,2]
    start=data_pos[0,0]+fix(s_pos*n_pix)
endif
;length=file_lines(file_spread)
;openr,1,file_spread
;data_spread=dblarr(5,length)
;readf,1,data_spread
;close,1
;x_spread=data_spread[3,*]
;xrms=2*data_spread[2,*]
;transalation in pixels to center of the detector default 257
im=make_array(nx,ny, /DOUBLE, value=0.0)

for x=0L,nx-1 do begin
  y=slope*x+intercept+y_shift
  y=fix(y)
; guassian spread: use when necessary 
; im[x,y+start]=ccd_count(x)
; npix=fix(x_spread[x])
; width=fix(xrms[x])
; xspread=psf_Gaussian( NPIXEL=npix, FWHM=width , /DOUBLE, /NORMALIZE);spread in cross dispesion direction guassian
; for j=0,npix-1 do begin
;   im[x,(y+start+j-(npix/2))]=xspread[j]*ccd_count(x)
; endfor

; Zemax value spread 
  ;disp_file='extra\xdisp_inter.txt'
  disp_file=file_spread
;use below file for the guassian functions
;disp_file='extra\xdisp_fnc.txt'

  len=file_lines(disp_file)
  openr,1,disp_file
  x_spread=dblarr(12,len)
  readf,1,x_spread
  close,1

; The spectrum is divided into 11 wavelength segments with varying cross-dispersion spread 
  if x le 180 then spread = x_spread[1,*]
  if x ge 180 and x le 360 then spread   = x_spread[2,*]
  if x ge 360 and x le 540 then spread   = x_spread[3,*]
  if x ge 540 and x le 720 then spread   = x_spread[4,*]
  if x ge 720 and x le 900 then spread   = x_spread[5,*]
  if x ge 900 and x le 1080 then spread  = x_spread[6,*]
  if x ge 1080 and x le 1260 then spread = x_spread[7,*]
  if x ge 1260 and x le 1440 then spread = x_spread[8,*]
  if x ge 1440 and x le 1620 then spread = x_spread[9,*]
  if x ge 1620 and x le 1800 then spread = x_spread[10,*]
  if x ge 1800 then spread               = x_spread[11,*]

  for j=0,len-1 do begin
    out=spread[j]*ccd_count(x)
    ;spread_value=out
    spread_value=out
    if ((y+start+x_spread[0,j]) lt ny-1 and (y+start+x_spread[0,j]) gt 0) then im[x,(y+start+x_spread[0,j])]=spread_value
  endfor

endfor
;stop
return,im
end