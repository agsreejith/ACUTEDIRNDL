; NAME:
;      CUTE_BACKGROUND
;
; PURPOSE:
;      Returns a 2D array that is the background frame for CUTE spectrograph;
; CALLING SEQUENCE:
;      Result = CUTE_BACKGROUND(nx, ny, token, ra, dec, jd, file_bg, file_sol, file_zod, wavelength )
;
; INPUTS:
;      nx,ny = The number of pixels in x and y for CUTE spectrograph detector
;      token = The values specifying if default background is to be loaded or not, y for default, n for calculating
;      Ra = The RA values of star, set to zero if default mode is used;
;      Dec = The Dec value of the star, set to zero if default mode is used
;      JD = The Julian date asosiated with the observation, set to zero if default mode is used
;      file_bg = Location of default background values
;      flie_sol = Location of solar spectrum
;      file_zod = Location of Zodiacal data
;      wavelength = wavelength array of cute
;      
;      
; OUTPUT:
;      Returns the CUTE background frame
; REQUIRES:
;      The background files stored at exact locations
;
; PROCEDURE:
;      Calculate the background based on zodiacal data or STIS high background 
;
; MODIFICATION HISTORY:
;      created 01.12.2017 by A. G. Sreejith

function cute_background, nx=nx, ny=ny, token=token, ra=ra, dec=dec, jd=jd, file_bg=file_bg, file_sol=file_sol, file_zod=file_zod, wavelength=wavelength

if not keyword_defined(nx       ) then nx        = 2048
if not keyword_defined(ny       ) then ny        = 515
if not keyword_defined(token    ) then token     = 'default'
if not keyword_defined(ra       ) then ra        = 0
if not keyword_defined(dec      ) then dec       = 0
if not keyword_defined(jd       ) then jd        = 0
if not keyword_defined(file_bg  ) then file_bg   = 'extra\background.txt'; in ergs/s/cm2/A/arcsec2 vs wavelength in A
if not keyword_defined(file_sol ) then file_sol  = 'extra\zod_spectrum.txt'
if not keyword_defined(file_zod ) then file_zod  = 'extra\zod_dist.txt'

background=dblarr(nx)

  if (token eq 'default') then begin
    length=file_lines(file_bg)
    data_bg=dblarr(2,length)
    openr,1,file_bg
    readf,1,data_bg
    close,1
    background_in=data_bg[1,*]*25*5.03e7*data_bg[0,*]; converting to photons/s/cm2/per pixel
    background=interpol(background_in,data_bg[0,*],wavelength,/SPLINE)
  endif else if(token eq 'calc') then begin
    ;background function
    ;background will depend on look direction and date
    ;Zodiacal Light in units of 10-8 W cm-2 sr-1 micron-1 at 5000 A. To convert into ph m-2 s-1 sr-1 A-1 at 5000 A, multiply by 252
    sunpos,jd,ra_s,dec_s
    euler,ra_s,dec_s,elb_s,ela_s,select=3
    length=file_lines(file_zod)
    data_zod=dblarr(12,length)
    openr,1,file_zod
    readf,1,data_zod
    close,1
    euler,ra,dec,elb,ela,select=3

    ela_h = fix(ela-ela_s)
    if (ela_h > 180) then ela_h = 360 - ela_h;
    elb_h = fix(elb)
    ;finding ph cm-2 s-1 sr-1 A-1 at 5000 A due to zodiacal
    i=0
    while(data_zod[i,0] lt ela_h) do i=i+1
    if (i gt 0) then begin
      if (( data_zod[i,0] - ela_h) gt  (ela_h - data_zod[i-1,0])) then i=i-1
    endif
    j=0
    while(data_zod[0,j] lt elb_h) do j++
    if (j gt 0) then begin
      if (( data_zod[0,j] - elb_h) gt  (elb_h - data_zod[0,j-1])) then j--
    endif
    zod_value=data_zod[i,j]

    ;The solar spectrum in photon units scaled such that the flux at 5000 A is equal to 10-8 W m-2 sr-1 micron-1.

    file_zods=file_sol
    length=file_lines(file_zods)
    data_zod=dblarr(2,length)
    openr,1,file_zods
    readf,1,data_zod
    close,1

    zod_spectrum=zod_value*data_zod[1,*]

    background=interpol(zod_spectrum,data_zod[0,*],wavelength,/SPLINE)
    ; to photons/s/cm2/pixel (use: 1 steradian = 1 rad2 = 3282.8 deg2 = 4.25 x 10^10 arcsec2, where 4p steradians = sphere.)
    background=background*25/4.25e10
  endif else begin
    print,'Invalid background input: check parameter file'
  endelse
cute_background=background


  return, cute_background
end