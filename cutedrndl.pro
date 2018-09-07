;//A CUTE data simulator.
;*******************************
; Copyright 2017 
;*This program is free software: you can redistribute it and/or modify
;*it under the terms of the GNU General Public License as published by
;*the Free Software Foundation, either version 3 of the License, or
;*(at your option) any later version.
;*
;*This program is distributed in the hope that it will be useful,
;*but WITHOUT ANY WARRANTY; without even the implied warranty of
;*MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;*GNU General Public License for more details.
;*
;*You should have received a copy of the GNU General Public License
;*along with this program.  If not, see <http://www.gnu.org/licenses/>.
; NAME:
;      CUTEDRNDL
;
; PURPOSE:
;      Creates CCD files for CUTE spectrograph;
; CALLING SEQUENCE:
;      CUTEDRNDL
;
; INPUTS:
;      parameter file in the current directory as the code, if else edit the code to provide the location of parameter file
;
; OUTPUT:
;      Returns raw, dark, bias, flat and background frames
; REQUIRES:
;      The requirements are provided in the parameter file and cute_readme.txt assosiated with the code.
;
; PROCEDURE:
;      Create RAW CCD frames fpr CUTE spectrograph based on temperature, radis and magnitude of a star;
; MODIFICATION HISTORY:
;      created 04.12.2017 by A. G. Sreejith
;      modified 20.12.2017  
;      modified 01.01,2018 
;      modofied 23.01.2018
;      modified 06.02.2018


FUNCTION RemoveRows, array, rows

  ; array -- A 2D array from which rows will be removed.
  ; rows -- A vector of row indices to remove from array.

  ; Need both positional parameters.
  IF N_Params() NE 2 THEN BEGIN
    Print, "Usage: 'RemoveRows, array, rowsToRemove'"
    RETURN, -1
  ENDIF

  ; The array must be 2D.
  ndims = Size(array, /N_DIMENSIONS)
  IF ndims NE 2 THEN BEGIN
    void = Dialog_Message('Array must be 2D.')
    Print, "Usage: 'RemoveRows, array, rowsToRemove'"
    RETURN, -1
  ENDIF

  ; The rows must be a vector.
  IF Size(rows, /N_DIMENSIONS) EQ 0 THEN rows = [rows]

  ; Find the dimensions of the array.
  dims = Size(array, /DIMENSIONS)

  ; Return the shortened array.
  RETURN, array[*, Where(~Histogram(rows, MIN=0, MAX=dims[1]-1), /NULL)]

END


pro cutedrndl

;for qwerty=4,29 do begin
;
;  #Configuration file for the data simulator
;  #provided here for easy reference
;  #Please note: comments and empty lines have to start with '#'
;  #Location of stellar input files
;  input_file      = 'D:\cutedrndl\'
;  #
;  #Stellar Parameters
;  #The parameters for the input star: Temperature in K, Radius in solar radii and V-Magnitude
;  #
;  stellar_temperature = 6100
;  stellar_radius    = 1.125
;  stellar_magnitude = 7.65
;  #
;  #Observing Parameters, RA and Dec in degrees, exposure time (exptime) in seconds, set JD to 0 for system time else specify JD
;  Ra          = 300.1792
;  Dec         = 22.7108
;  JD          = 0
;  exptime       = 300
;  #
;  #Instrument Parameters :
;  #Position of the star in slit : 3 options avaliable (0',±4',±8' from slit center), input 0,4 or 8 respectively
;  #Spectral resolution of the spectrograph in Å.
;  #
;  slit_position   = 0
;  spectral_resolution = 0.8
;  #
;  #Wavelength file location for CUTE
;  wave_file     = 'extra\wavelength.txt'
;  #
;  #Effective area file location for CUTE
;  file_eff      = 'extra\eff_area.txt'
;  #
;  #Spectrum parameter file location for CUTE
;  #Specifies the position on the detector, the slope and y intercept of the the spectrum.
;  #Default file has 3*3 array with each row coresponding to the coresponding position on the slit from center(0)
;  #
;  file_pos      = 'extra\position.txt'
;  #Specify the file contating information regarding spread of the spectrum in cross dispresion axis
;  #
;  file_spread     = 'extra\xdisp_inter.txt'
;  #
;  #Background information
;  #Select type of background: Options default or calc,
;  #default: STIS high background values. calc: Calculated from zodiacal light data
;  #Include stars, to include stars chose 1 else 0
;  type_bckg     = 'default'
;  stars       = 1
;  #
;  #file location for stellar parameters for different stellar type
;  stellar_params    = 'extra\stellar_param.txt'
;  #file locations for background calculation
;  bckg_default='extra\background.txt'
;  #
;  #file locations for zodiacal calculation: solar spectrum and zodiacal distribution
;  file_sol='extra\zod_spectrum.txt'
;  file_zod  = 'extra\zod_dist.txt'
;  #
;  #CCD parameters
;  #for 300sec, change values for coresponding exposure time
;  #Readout noise in e/pixel
;  readout_noise   = 12.25
;  bias_value      = 8
;  dark_value      = 3.6
;  ccd_gain      = 1.2
;  x_pixels      = 2048
;  y_pixels      = 515
;  #
;  #Output location
;  file_out      = 'output\'
;  #
;  #Spacecraft parameters
;  #To include effect of jitter in the simulator  chose 1 for jitter_sim and 0 to skip jitter
;  #Jitter of spacecraft RMS value, Rotation wrt to RA axis
;  jitter        = 7.2
;  rotation      = 0
;  jitter_sim      = 1
;  #
;  #Transit parameters
;  #Orbita period is in days, inclination in degrees, semi major axis and
;  #planet radius as a fraction of stellar radius
;  limb_file     = 'D:\cutedrndl\limb_darkening\phoenix\'
;  impact_parameter  = 0.3877
;  orbital_period    = 2.788
;  orbital_inclination = 87.18
;  rlt_semi_major_axis = 7.88
;  #Either provide planet radius as one number of the location of the file containing variation of
;  #planet radius (ratio of planet radius to stellar radius) with wavelength in CUTE band.
;  planet_radius   = 0.113
;  #planet_radius_file  = 'extra\radius.txt'

;Main code

close,/all
;logfile


;Reading from configuration file
infile=gm_read_textstructure("cutedrndl_parameters_hd209.txt");edit the configuration information above and comment this line to run without parameter file.
file=infile.input_file
t_star=double(infile.stellar_temperature)
r_star=double(infile.stellar_radius)
M_star=double(infile.stellar_magnitude)
s_pos=double(infile.slit_position)
if (tag_exist(infile,'spectral_resolution') eq 1) then fwhm=float(infile.spectral_resolution)
file_wave=infile.wave_file
file_eff=infile.file_eff
file_pos=infile.file_pos
file_spread=infile.file_spread
key=infile.type_bckg
file_bg=infile.bckg_default
file_sol=infile.file_sol
file_zod=infile.file_zod
file_out=infile.file_out
r_noise_in=double(infile.readout_noise)
exptime=long(infile.exptime)
r_time=long(infile.read_time)
b_value=double(infile.bias_value)
dark_value=double(infile.dark_value)
G=double(infile.ccd_gain)
ra=double(infile.Ra)
dec=double(infile.Dec)
jd=double(infile.JD)
xjitter=double(infile.jitter_in_x)
yjitter=double(infile.jitter_in_y)
zjitter=double(infile.jitter_in_z)
nx=long(infile.x_pixels)
ny=long(infile.y_pixels)
n_cosmic=long(infile.cosmic)
;read out noise in eloectrons
r_noise=r_noise_in/G ;use this if r_noise in in counts
r_noise=r_noise_in
s_wid=infile.slit_width
s_len=infile.slit_length
pl_scale=double(infile.plate_scale)

;file_out='output\hd209\bz\orbit'+String(qwerty, Format='(I02)')+'\'
;slit position check
if (abs(s_pos) gt double(s_len/120)) then begin
  print,'Slit position outside slit, Readjust the parameters and run again'
  goto, ending
endif
;get system time
if jd eq 0 then t=systime(/JULIAN, /UTC) else t = jd

;read input wavelength file
w_length=file_lines(file_wave)
wavelength=dblarr(w_length)
ccd_flux=dblarr(nx)
ccd_wave=dblarr(nx)
openr,1,file_wave
readf,1,wavelength
close,1

;if resolution determine the the wavelength and number of pixels in x
if (infile.developer_mode eq 1) then begin
  start_wave=double(infile.start_wave)
  wavelength=make_array(DIMENSION=nx, INCREMENT=fwhm/2, /INDEX, /NOZERO, START=start_wave) 
  openu,1,file_wave
  for i=0L,n_elements(wavelength)-1 do printf,1, wavelength[i]
  close,1
endif

print,t
t1=systime(/SECONDS)
;convert stellar parameters to photons at earth
photons=cute_photons(infile) 
photons_star=photons[1,*]
wave=photons[0,*]
photons_star=reform(photons_star,n_elements(photons_star))
wave=reform(wave,n_elements(wave))
writecol,file_out+'photon_input.txt',wave,photons_star,fmt='(2(F17.9,1x))
;stop
;jitter test, creates a gaussian of specified sigma value to check the effect of jitter 
;jitter_flux=make_array(n_elements(photons_star),value=1)
;gauss_flux=gaussian(wave,[1000,2900,0.0213])
;st=where(wave eq 2850)
;en=where(wave eq 2950)
;jitter_flux=jitter_flux+gauss_flux
;wave_print=wave[st:en]
;jitter_print=jitter_flux[st:en]
;photons_star=jitter_flux
;writecol,file_out+'jitter_input.txt',wave_print,jitter_print,fmt='(2(F17.9,1x))


st=value_locate(wave, 2000)
en=value_locate(wave, 7000)
wave_new=wave[st:en]
photons_star_new=photons_star[st:en]
;for nonlinear dispersion
if (tag_exist(infile,'dispersion_file') eq 1) then begin
  ;if the spectrograph has non linear dispersion then
  dis_len=file_lines(infile.dispesrsion_file);contains wavelength and FWHM info
  openr,1,infile.dispesrsion_file
  readf,1,dispersion_data
  close,1
  wave_disp=dispersion_data[0,*]
  fwhm_disp=dispersion_data[1,*]
  fwhm=interpol(fwhm_disp,wave_disp,wave_new)
  ;wavelength resolution of CUTE
  if (s_pos ge -2 && s_pos le 2) then fwhm=fwhm
  if (s_pos gt 2 && s_pos le 8) then fwhm=fwhm*3/2
  if (s_pos gt 8) then fwhm=fwhm*2
  if (s_pos lt -2 && s_pos ge -8) then fwhm=fwhm*3/2
  if s_pos lt -8 then fwhm=fwhm*2
 
  hwhm=fwhm/2
  f_flux=dblarr(n_elements(wave_new))
  for i=0,n_elements(wave_new)-1 do begin
    wave_s=dblarr(n_elements(wave_new))
    flux_s=dblarr(n_elements(photons_star_new))
    wave_s=wave_new
    flux_s[i]=photons_star_new[i]
    gauss_flux=gaussbroad(wave_s,flux_s,hwhm[i])
    f_flux+=gauss_flux
  endfor
  smoothedflux=f_flux  
 endif else begin ;linear dispersion value set by spectral resolution
  ;wavelength resolution of CUTE
  if (s_pos ge -2 && s_pos le 2) then fwhm=fwhm
  if (s_pos gt 2 && s_pos le 8) then fwhm=fwhm*3/2
  if (s_pos gt 8) then fwhm=fwhm*2
  if (s_pos lt -2 && s_pos ge -8) then fwhm=fwhm*3/2
  if s_pos lt -8 then fwhm=fwhm*2 
  hwhm=fwhm/2

st1=value_locate(wave_new, 2518.50)
en1=value_locate(wave_new, 3331.10)

  ;convolution with instrument response
  smoothedflux=gaussbroad(wave_new,photons_star_new,hwhm)
  ;smoothedflux = gaussfold(wave, photons_star, fwhm, LAMMIN=2518.50, LAMMAX=3331.10)
endelse
wave_res=dblarr(nx/2)
for i=0,nx-2,2 do wave_res[i/2]=mean(wavelength[i:i+1])

;interpolate and trim to detector size
;ccd_flux=interpol(smoothedflux,wave_new,wavelength,/SPLINE) ;for jitter
ccd_flux=interpol(smoothedflux,wave_new,wave_res,/SPLINE); linear interpolation
ccd_wave=wavelength

;read effective area and interpolate to cute wavelength
length=file_lines(file_eff)
eff_area=dblarr(5,length)
openr,1,file_eff
readf,1,eff_area
close,1
aeff=interpol(eff_area[4,*],eff_area[0,*],wave_res,/SPLINE); linear interpolation

ccd_count1=dblarr(nx/2) 
;Effective area/QE
ccd_count1=ccd_flux;*aeff*fwhm ;comment *aeff*fwhm during jitter test
plot,wave_res,ccd_count1
writecol,file_out+'cute_counts.txt',wave_res,ccd_count1*300.
;write_png,file_out+'cute_count.png',TVRD(/TRUE)
ccd_count=dblarr(nx)

for i=0,nx-1 do begin
 ccd_count[i]=ccd_count1[i/2]/2
 ;print,ccd_count1[i/2],ccd_count[i]
endfor
;for jitter
;ccd_count=ccd_flux
;writecol,file_out+'jitter_cute.txt',wavelength,ccd_count,fmt='(2(F17.9,1x))
writecol,file_out+'cute_ccd_counts.txt',wavelength,ccd_count*300.
;stop
;oplot,wave_res,ccd_count,color=cgcolor('red')
;stop
; y shift in CCD coresponding to different slit positions
;yshift=s_pos*60*2.5;input position is in arcmin converting it to pixel sclae 2.5" is the plate sclae of CUTE 
yshift=0;shift in position is now taken care by cute_specspread
;if s_pos eq 0 then yshift=0
;if s_pos eq 4 then yshift=96
;if s_pos eq 8 then yshift=192
;if s_pos eq -4 then yshift=-96
;if s_pos eq -8 then yshift=-192

;variables for total bias, dark, flat etc
t_bias=dblarr(nx,ny)
t_dark=dblarr(nx,ny)
t_flat=dblarr(nx,ny)
;defining cosmic ray frame
ray=make_array(nx,ny,value=0.0,/DOUBLE)

;bias
;generate bias frames
for k=0,4 do begin
  rand=fix(20*randomu(seed,1))
  ran_val=rand[0]
  for rd=0,ran_val do bias1=r_noise*randomn(seed,nx*ny*2)+b_value
  new_r=randomu(seed,1)
  if new_r gt 0.5 then new_bias=bias1[0:nx*ny-1]
  if new_r le 0.5 then new_bias=bias1[nx*ny:(nx*ny*2)-1]
  bias= reform(new_bias,[nx,ny])
  ;bias=cute_random(nx,ny,r_noise,b_value)
  ;stop
  sxaddpar, bihdr, 'Time_in_JD', t
  sxaddpar, bihdr, 'Readout_noise', r_noise
  sxaddpar, bihdr, 'Bias_value', b_value
  sxaddpar, bihdr, 'Mean', mean(bias)
  sxaddpar, bihdr, 'Median ', median(bias)
  sxaddpar, bihdr, 'Max_value', max(bias)
  sxaddpar, bihdr, 'Min_Value', min(bias)
  bias=bias*G
  ;write the fits file
  mwrfits,bias,file_out+'bias'+String(k, Format='(I05)') +'.fits',bihdr,/create
  t_bias=t_bias+bias
endfor
bias_num=bias
;average bias
bias=(t_bias/5)


for k=0,4 do begin
  ;dark
  rand=fix(20*randomu(seed,1))
  ran_val=rand[0]
  for rd=0,ran_val do bias1=r_noise*randomn(seed,nx*ny*2)+b_value
  new_r=randomu(seed,1)
  if new_r gt 0.5 then new_bias=bias1[0:nx*ny-1]
  if new_r le 0.5 then new_bias=bias1[nx*ny:(nx*ny*2)-1]
  bias= reform(new_bias,[nx,ny])
  ;bias=cute_random(nx,ny,r_noise,b_value)
  d=dblarr(nx*ny*2)
  rand=fix(20*randomu(seed,1))
  ran_val=rand[0]
  for rd=0,ran_val do d=0.001*randomn(seed,nx*ny*2)
  new_r=randomu(seed,1)
  if new_r gt 0.5 then new_dark=d[0:nx*ny-1]
  if new_r le 0.5 then new_dark=d[nx*ny:(nx*ny*2)-1]
  new_dark=new_dark+dark_value; dark of 3.6 e/pixel
  dark= reform(new_dark,[nx,ny]); in e/pixel
   ;dark=cute_random(nx,ny,0.001,dark_value)
   ;stop
  if n_cosmic ne 0 then ray=cosmic(n_cosmic,72000,nx,ny)
  dark=(dark*exptime)+bias
  dark2=dark+ray
  dark=dark*G
  dark2=dark2*G
  ;writing dark header
  sxaddpar, dhdr, 'Time_in_JD', t
  sxaddpar, dhdr, 'Input_dark_value', dark_value
  sxaddpar, dhdr, 'Mean', mean(dark)
  sxaddpar, dhdr, 'Median ', median(dark)
  sxaddpar, dhdr, 'Max_value', max(dark)
  sxaddpar, dhdr, 'Min_Value', min(dark)
  sxaddpar, dhdr, 'Exposure_Time_in_seconds', exptime
  
  ;write the fits file
  mwrfits,dark,file_out+'dark'+String(k, Format='(I05)') +'.fits',dhdr,/create
  mwrfits,dark2,file_out+'dark_wcosmic'+String(k, Format='(I05)') +'.fits',dhdr,/create
  t_dark=t_dark+dark
endfor
dark_n=dark
;average dark
dark=(t_dark/5.)

;creating background frame
backg=cute_background(wavelength=wavelength,token=key, file_bg=file_bg, file_sol=file_sol, file_zod=file_zod)

if(key eq 'default') then token='Default: From STIS backgrounds'
if(key eq'calc') then token= 'Calculated from Zodiacla data'
aeff2=interpol(eff_area[4,*],eff_area[0,*],wavelength,/SPLINE);
scaled_spectrum=backg*aeff2

;adding background to create an image array
background=dblarr(nx,ny)

for k=0,ny-1 do begin
  background(*,k)=scaled_spectrum
endfor
background=background*exptime
;background header
sxaddpar, bhdr, 'Time_in_JD', t
sxaddpar, bhdr, 'RA_of_target', ra
sxaddpar, bhdr, 'Dec_of_target', dec
sxaddpar, bhdr, 'Julian_date', jd
sxaddpar, bhdr, 'Type_of_background', token
sxaddpar, bhdr, 'Spectral_resolution', fwhm
sxaddpar, bhdr, 'Exposure_Time_in_seconds', exptime
;write the fits file
mwrfits,background,file_out+'background.fits',bhdr,/create

; checking for background stars and calculating background stars
im_bs=dblarr(nx,ny)
if infile.stars eq 1 then im_bs=background_star(infile)

;adding scatter
im_sc=dblarr(nx,ny)
if infile.scatter eq 1 then im_sc=cute_scatter(nx,ny,tele_area,scatter,reflection)
im_sc=im_sc*exptime

;generating flats
for k=0,4 do begin
  ;flat
  ;generate flatframe
  ;creating aflat frame
  rand=fix(20*randomu(seed,1))
  ran_val=rand[0]
  for rd=0,ran_val do bias=r_noise*randomn(seed,nx*ny*2)+b_value
  new_r=randomu(seed,1)
  if new_r gt 0.5 then new_bias=bias1[0:nx*ny-1]
  if new_r le 0.5 then new_bias=bias1[nx*ny:(nx*ny*2)-1]
  bias= reform(new_bias,[nx,ny])
  ;bias=cute_random(nx,ny,r_noise,b_value)
  d=dblarr(nx*ny*2)
  rand=fix(20*randomu(seed,1))
  ran_val=rand[0]
  for rd=0,ran_val do d=0.001*randomn(seed,nx*ny*2)
  new_r=randomu(seed,1)
  if new_r gt 0.5 then new_dark=d[0:nx*ny-1]
  if new_r le 0.5 then new_dark=d[nx*ny:(nx*ny*2)-1]
  new_dark=new_dark+dark_value; dark of 3.6 e/s/pixel
  dark= reform(new_dark,[nx,ny]); in e/pixel
  ;dark=cute_random(nx,ny,0.001,dark_value)
  rand=fix(20*randomu(seed,1))
  ran_val=rand[0]
  for rd=0,ran_val do flat1=5*randomn(seed,nx*ny*2)+50000L
  new_r=randomu(seed,1)
  if new_r gt 0.5 then new_flat=flat1[0:nx*ny-1]
  if new_r le 0.5 then new_flat=flat1[nx*ny:(nx*ny*2)-1]
  flat=reform(new_flat,[nx,ny])
  ;flat=cute_random(nx,ny,5,50000L)
  if n_cosmic ne 0 then ray=cosmic(n_cosmic,72000,nx,ny)
  flat_n=flat/50000L
  flat=flat+dark+bias
  flat2=flat+ray
  flat=flat*G
  flat2=flat2*G
  ;write the fits file
  sxaddpar, fhdr, 'Time_in_JD', t
  sxaddpar, fhdr, 'Mean', mean(flat)
  sxaddpar, fhdr, 'Median ', median(flat)
  sxaddpar, fhdr, 'Max_value', max(flat)
  sxaddpar, fhdr, 'Min_Value', min(flat)
  mwrfits,flat,file_out+'flat'+String(k, Format='(I05)') +'.fits',fhdr,/create
  mwrfits,flat2,file_out+'flat_wcosmic'+String(k, Format='(I05)') +'.fits',fhdr,/create
  t_flat=t_flat+flat
endfor  
;average flat
flat=(t_flat/5.)

im=make_array(nx,ny, /DOUBLE, value=0.0)
;calculating the scaling value based on transit parameters
if (tag_exist(infile,'transit') eq 0)then begin 
  transit_lightcurve,infile,scale,new_time,mid_time 
endif else begin  ;reading the planet transmission spectrum from file 
  file=infile.transit
  transit_len=file_lines(infile.transit)
  ; Determine the number of rows in the file.
  rows = file_lines(file)
  ; Determine the number of colums in the file by reading
  ; the first line and parsing it into column units.
  openr, lun, file, /Get_Lun
  line = ""
  readf, lun, line
  ; Find the number of columns in the line.
  cols = n_elements(StrSplit(line, /RegEx, /Extract))
  stop
  if cols eq 2 then begin ;if transit light curve is constant for all wavelengths
    data=dblarr(2,transit_len)
    openr,11,infile.transit
    readf,11,data
    close,11
    in_time=data[0,*];asuming time in seconds
    scale1=data[1,*]
    ;scale=interpol(sacel1,in_wave,ccd_wave)
    obs_time=in_time[n_elements(in_time)-1]-in_time[1]
    rounding = obs_time mod (exptime+r_time)
    reminder=exptime+r_time-rounding
    end_time = obs_time-rounding
    steps=value_locate(in_time, exptime+r_time)
    en= value_locate(in_time,end_time)
    j=0
    exp_step=value_locate(in_time,exptime)
    new_time=dblarr((en/steps)-1)
    scale=dblarr(nx,(en/steps)-1)
    jj=0
      for j=0, n_elements(scale[0,*])-1 do begin
        ;print,jj,j,jj+exp_step
        new_impct=scale1[jj:jj+exp_step]
        scale[*,j]=mean(new_impct)
        jj=jj+steps
        new_time[j]=in_time[jj]
        loc=where(scale_value eq min(scale_value))
        mid_time=in_time[loc]
      endfor
  endif else begin
    ; Create a variable to hold the data.
    data = dblarr(cols, rows)
    ; Rewind the data file to its start.
    point_lun, lun, 0
    ; Read the data.
    readf, lun, data
    free_lun, lun
    in_wave=data[*,0]
    in_time=data[0,*];asuming time in seconds
    scale_value=dblarr(cols-1,rows-1)
    scaling=dblarr(nx,rows-1)
    for m=1, cols-1 do scale_value[m-1,*]=data[m,*]
    for mn=1,rows-1 do scaling[*,mn-1]=interpol(scale_value[*,mn-1],in_wave,wavelength)
    obs_time=in_time[n_elements(in_time)-1]-in_time[1]
    rounding = obs_time mod (exptime+r_time)
    reminder=exptime+r_time-rounding
    end_time = obs_time-rounding
    steps=value_locate(in_time, exptime+r_time)
    en= value_locate(in_time,end_time)
    j=0
    exp_step=value_locate(in_time,exptime)
    new_time=dblarr((en/steps)-1)
    scale=dblarr(nx,(en/steps)-1)
    for k=0,nx-1 do begin
        jj=0
        for j=0, n_elements(scale[0,*])-1 do begin
           ;print,jj,j,jj+exp_step
            new_impct=scaling[k,jj:jj+exp_step]
            scale[k,j]=mean(new_impct)
            jj=jj+steps
            new_time[j]=in_time[jj]
            loc=where(scale_value eq min(scale_value))
            mid_time=in_time[loc]  
        endfor
    endfor
  ;scale=interpol(sacel1,in_wave,ccd_wave) 
  endelse
endelse  
m_t=mid_time/86400.0
Tmid=t+m_t[0]

if (tag_exist(infile,'transit_mid_time') eq 1)then Tmid=double(mid_time)  
;define header of image
writecol,file_out+'input_light curve.txt',new_time,scale(0,*),fmt='(2(F17.7,1x))
t0=t
if (tag_exist(infile,'transit_mid_time') eq 1)then t0=new_time[0] 

;implementing occultations
if infile.sc_orbit_period eq 0 then new_scale=scale else begin 
  ;cgplot,new_time,scale(0,*)*exptime,yrange=[240,310]
  time_factor=1.
  time_fact=0.
  orbit_period=double(infile.sc_orbit_period)*60.
  expos_fact=make_array((n_elements(new_time)),value=1.0)
  occ_time=double(infile.ocl_time)
  time_obs=exptime+r_time
  gap_dur=occ_time*60.
  orb_length=orbit_period-gap_dur
  rand=fix(20*randomu(seed,1))
  ran_val=rand[0]
  for rd=0,ran_val do gap=fix(time_obs+orb_length*RANDOMU(seed,1))
  gap_loc1=gap
  time_orbit=(new_time-t0)*86400.
  orbit_dur=time_orbit[n_elements(time_orbit)-1]
  n_orbits=1+fix(orbit_dur/orbit_period)
  gap_loc=dblarr(n_orbits)
  near=dblarr(n_orbits)
  index=intarr(n_orbits)
  end_loc=intarr(n_orbits)
  for i=0,n_orbits-1 do begin
    gap_loc[i]= gap_loc1+orbit_period*i
    near[i] = Min(Abs(time_orbit - gap_loc[i]), ind)
    diff = Min(Abs(time_orbit - (gap_loc[i]+gap_dur)), endloc1)
    index[i]=ind
    end_loc[i]=endloc1
  endfor
;if time_orbit gt orbit_period then time_orbit=time_orbit-orbit_period*fix(time_orbit/5400.)
  point=dblarr(n_orbits)
  time_fact=dblarr(n_orbits)
  for i=0,n_orbits-1 do begin
    if time_orbit[index[i]] gt gap_loc[i] then point[i]=time_orbit[index[i]-1] else point[i]=time_orbit[index[i]-1]
    if gap_loc[i]-time_orbit[index[i]] lt 0 then time_fact[i]= ((time_orbit[index[i]]-gap_loc[i]))/exptime
  endfor
  for i=0, n_orbits-1 do expos_fact[index[i]-1]= expos_fact[index[i]-1]-time_fact[i]
  new_length=(n_elements(time_orbit)-total(end_loc-index))
  removables=intarr(1)
  upd_time=dblarr(new_length)
  upd_scale=dblarr(nx,new_length)
  upd_exptime=dblarr(new_length)
  for i=0,n_orbits-1 do begin
    size=end_loc[i]-index[i]+1
    array=make_array(size,start=index[i],/INDEX)
    removables= [removables,array]
  endfor
  remove,0,removables
  remove,removables,new_time
  remove,removables,expos_fact
  new_scale=RemoveRows(scale,removables)
  ;cgoplot,new_time,new_scale(0,*)*exptime,psym=2,color='red'
  ;cgoplot,new_time,new_scale(0,*)*exptime*expos_fact,psym=2,color='blue'
  ;write_png,file_out+'occl.png',TVRD(/TRUE)
endelse
;stop
for k=0, n_elements(new_scale[0,*])-1 do begin ;
;  inc=k*exptime/86400.0
;  t=t0+inc
  if (tag_exist(infile,'transit_mid_time') eq 1)then t=new_time[k] 
  scaled_ccd_count=new_scale[*,k]*ccd_count
  ;cute_orbit,
  
  ;if t-t0 
  ;defien position on detector
  im=cute_specspread(file_pos,file_spread,scaled_ccd_count,nx,ny,s_pos,0,yshift,pl_scale)
  ;mwrfits,im,file_out+'1sec_image_wtc'+String(k, Format='(I05)') +'.fits',hdr1,/create
  im_wtc=im*exptime*G
  ;adding backgroubd stars
  im_t=im+im_bs
  delta_x=xjitter
  delta_y=yjitter
  delta_z=zjitter
  ;exptime=fix(exptime*expos_fact[k])
  if infile.jitter_sim eq 1 then im_j=jitter(im_t,delta_x,delta_y,delta_z,exptime,pl_scale) else if infile.jitter_sim eq 0 then im_j=im*exptime
  im_jl=im_j ;if no linearity
  ;adding the effect of linearity
  image_linearity,im_j,im_jl
  ;stop
  if infile.systematics eq 1 then im_jls=cute_systematics(infile,t,t0,im_jl) else im_jls=im_jl ;adding cute systematics
  
  sxaddpar, hdr1, 'Time_in_JD', t
  sxaddpar, hdr1, 'Star_Tem', t_star
  sxaddpar, hdr1, 'Star_Rad', r_star
  sxaddpar, hdr1, 'Magnitude ', M_star
  sxaddpar, hdr1, 'Slit_Position', s_pos
  sxaddpar, hdr1, 'Spectral_resolution', fwhm
  sxaddpar, hdr1, 'Exp_Time', exptime
  ;write the fits file
  sxaddpar, hdr, 'obs_time', t,FORMAT='F17.7'
  sxaddpar, hdr, 'Star_Temp', t_star
  sxaddpar, hdr, 'Star_Rad', r_star
  sxaddpar, hdr, 'Magnitude', M_star
  sxaddpar, hdr, 'Slit_Pos', s_pos
  sxaddpar, hdr, 'Spectral', fwhm
  sxaddpar, hdr, 'Exp_Time', exptime
  sxaddpar, hdr, 'mid_time', Tmid[0],FORMAT='F17.7'
 

  
 ;mwrfits,im_wtc,file_out+'image_wtc'+String(k, Format='(I05)') +'.fits',hdr1,/create
;bias, flat and dark check
;  im_wtc_num=im_wtc+bias_num 
;  mwrfits,im_wtc_num,file_out+'image_wtc_bias'+String(k, Format='(I05)') +'.fits',hdr1,/create
;  im_wtc_drk=im_wtc+dark_n
;  im_wtc_flt=im_wtc/flat_n
;  mwrfits,im_wtc_drk,file_out+'image_wtc_drk'+String(k, Format='(I05)') +'.fits',hdr1,/create
;  mwrfits,im_wtc_flt,file_out+'image_wtc_flt'+String(k, Format='(I05)') +'.fits',hdr1,/create
;  stop
  
  im_n =dblarr(nx,ny)
  ;readoutnoise
  noi=dblarr(nx,ny)
  ;add Photon noise
  im_n= im_jls+sqrt(im_jls)
  noi=sqrt(im_jls)
    
  im_nb=im_n+background+im_sc
  ;defining flat
  rand=fix(20*randomu(seed,1))
  ran_val=rand[0]
  for rd=0,ran_val do flat1=0.0001*randomn(seed,nx*ny*2)+1
  new_r=randomu(seed,1)
  if new_r gt 0.5 then new_flat=flat1[0:nx*ny-1]
  if new_r le 0.5 then new_flat=flat1[nx*ny:(nx*ny*2)-1]
  flat=reform(new_flat,[nx,ny])
  ;flat=cute_random(nx,ny,0.0001,1)
  
  ;im_nbf=im_nb/flat ;redundant now
  ;defining bias
  rand=fix(20*randomu(seed,1))
  ran_val=rand[0]
  for rd=0,ran_val do bias1=r_noise*randomn(seed,nx*ny*2)+b_value
  new_r=randomu(seed,1)
  if new_r gt 0.5 then new_bias=bias1[0:nx*ny-1]
  if new_r le 0.5 then new_bias=bias1[nx*ny:(nx*ny*2)-1]
  bias= reform(new_bias,[nx,ny])
  ;bias=cute_random(nx,ny,r_noise,b_value)
  ;im_nbfb=im_nbf+bias
  
  ;defining dark
  d=dblarr(nx*ny*2)
  rand=fix(20*randomu(seed,1))
  ran_val=rand[0]
  for rd=0,ran_val do d=0.001*randomn(seed,nx*ny*2)
  new_r=randomu(seed,1)
  if new_r gt 0.5 then new_dark=d[0:nx*ny-1]
  if new_r le 0.5 then new_dark=d[nx*ny:(nx*ny*2)-1]
  new_dark=new_dark+dark_value; dark of 3.6 e/pixel
  dark= reform(new_dark,[nx,ny])
  ;dark=cute_random(nx,ny,0.001,dark_value)
  im_nbfbd=(im_nb+(dark*exptime)+bias)*flat
  ;;to counts where G is the gain
  im_final=im_nbfbd*G
  
  ;cosmic ray generation
  ;print,k
  if n_cosmic ne 0 then ray=cosmic(n_cosmic,72000,nx,ny)
  im_final=im_final+ray
  ;mwrfits,ray,file_out+'cosmic'+String(k, Format='(I05)') +'.fits',hdr,/create
  ;write the fits file
  mwrfits,im_final,file_out+'raw_image'+String(k, Format='(I05)') +'.fits',hdr,/create
  inc=double((exptime+r_time)/86400.0)
  t=t+inc

endfor
im_bs=im_bs*G*exptime
  mwrfits,im_bs,file_out+'background_star.fits',hdr,/create
t2=systime(/SECONDS)  
  print,'Exitng the simulator'
  openw,12,'cutedrndl_log_at_'+String(k, Format='(I05)') +'.txt'
  printf,12,'Input parameters'
  printf,12,'Stellar Temperature:',infile.stellar_temperature
  printf,12,'Stellar Radius:',infile.stellar_radius
  printf,12,'Stellar Magnitude:',infile.stellar_magnitude
  printf,12,'Position in slit:',infile.slit_position
  printf,12,'Exposure time:',exptime
  printf,12,'Background calculated from',token
  if infile.stars eq 1 then printf,12,'Background stars included from Vizier'
  if (tag_exist(infile,'transit_mid_time') eq 1)then printf,12,'Transit mid time set from user input' else printf,12,'Transit mid time calculated from JD'
t2=systime(/SECONDS)
ending:
print,'Time taken by simulator in seconds:',t2-t1
;udefine all variable
;delvar,bias,dark,flat,im,im_j,im_jl,im_jln,im_n,im_nb,im_nbfbd,im_final,t,k,st,en,ccd_count
;endfor
end