; NAME:
;      CUTE_POSITION
;
; PURPOSE:
;      Returns counts in photons/s/cm2 for a star of certain magnitude, radius and temperature;
; CALLING SEQUENCE:
;      Result = CUTE_PHOTONS(file, t_star, r_star, m_star)
;
; INPUTS:
;      file = The location of stellar spectra.
;      t_star = The temperature of the star, 4100 if keyword is absent.
;      r_star = The radius of  the star (in solar radius),set to one solar radius if keyword is absent;
;      m_star = The visula magnitude of the star, set to 10 is keyword is absent.
;
; OUTPUT:
;      Returns counts in photons/s/cm2
; REQUIRES:
;      The input stellar spectrum files stored at exact locations
;
; PROCEDURE:
;      Calculate the counts at earth in photons/s/cm2 based on temperature, radis and magnitude of a star;
; MODIFICATION HISTORY:
;      created 01.12.2017 by A. G. Sreejith
;      extinction added xx.12.2017 by A. G. Sreejith

function cute_position, file_pos, file_spread,

  if not keyword_defined(file_pos   ) then file_pos    = 'extra\position.txt'
  if not keyword_defined(file_spread) then file_spread ='extra\cuteback.txt'


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
openr,1,file_spread
data_spread=dblarr(5,2048)
readf,1,data_spread
close,1
x_spread=data_spread[3,*]
xrms=2*data_spread[2,*]
;transalation in pixels to center of the detector default 257
im=dblarr(2048,515)

for x=0,2047 do begin
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

end