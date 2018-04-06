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


pro cutedrndl

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
infile=gm_read_textstructure("cutedrndl_parameters.txt");edit the configuration information above and comment this line to run without parameter file.
file=infile.input_file
t_star=double(infile.stellar_temperature)
r_star=double(infile.stellar_radius)
M_star=double(infile.stellar_magnitude)
s_pos=double(infile.slit_position)
fwhm=float(infile.spectral_resolution)
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
;slit position check
if (abs(s_pos) gt double(s_len/120)) then begin
  print,'Slit position outside slit, Redajust the parameters and run again'
  goto, ending
endif
;get system time
if jd eq 0 then t=systime(/JULIAN, /UTC) else t = jd

print,t
t1=systime(/SECONDS)
;convert stellar parameters to photons at earth
photons=cute_photons(infile) 
;printf,12,'Photons calculated' 
photons_star=photons[1,*]
wave=photons[0,*]
photons_star=reform(photons_star,n_elements(photons_star))
wave=reform(wave,n_elements(wave))

;jitter test, creates a gaussian of specified sigma value to check the effect of jitter 
;jitter_flux=make_array(n_elements(photons_star),value=1)
;gauss_flux=gaussian(wave,[1000,2900,1.8])
;st=where(wave eq 2850)
;en=where(wave eq 2950)
;jitter_flux=jitter_flux+gauss_flux
;wave_print=wave[st:en]
;jitter_print=jitter_flux[st:en]
;photons_star=jitter_flux
;writecol,file_out+'jitter_input.txt',wave_print,jitter_print

;wavelength resolution of CUTE
if (s_pos ge -2 && s_pos le 2) then fwhm=fwhm
if (s_pos gt 2 && s_pos le 8) then fwhm=fwhm*3/2
if (s_pos gt 8) then fwhm=fwhm*2
if (s_pos lt -2 && s_pos ge -8) then fwhm=fwhm*3/2
if s_pos lt -8 then fwhm=fwhm*2 

st=value_locate(wave, 2000)
en=value_locate(wave, 3500)
wave_new=wave[st:en]
photons_star_new=photons_star[st:en]
hwhm=fwhm/2

;convolution with instrument response
smoothedflux=gaussbroad(wave_new,photons_star_new,hwhm)
;smoothedflux = gaussfold(wave, photons_star, fwhm, LAMMIN=2518.50, LAMMAX=3331.10)


w_length=file_lines(file_wave)
wavelength=dblarr(w_length)
ccd_flux=dblarr(nx)
ccd_wave=dblarr(nx)
openr,1,file_wave
readf,1,wavelength
close,1

;interpolate and trim to detector size
ccd_flux=interpol(smoothedflux,wave_new,wavelength,/SPLINE); linear interpolation
ccd_wave=wavelength


;read effective area and interpolate to cute wavelength
length=file_lines(file_eff)
eff_area=dblarr(11,length)
openr,1,file_eff
readf,1,eff_area
close,1
aeff=interpol(eff_area[4,*],eff_area[0,*],ccd_wave,/SPLINE); linear interpolation

ccd_count=dblarr(nx) 
;Effective area/QE
ccd_count=ccd_flux*aeff

; y shift in CCD coresponding to different slit positions
;yshift=s_pos*60*2.5;input position is in arcmin converting it to pixel sclae 2.5" is the plate sclae of CUTE 
yshift=0;shift in position is taken care by cute_specspread
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
  bias=r_noise*randomn(seed,nx*ny)+b_value
  bias= reform(bias,[nx,ny])
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
  bias=r_noise*randomn(seed,nx*ny)+b_value
  bias= reform(bias,[nx,ny])
  d=dblarr(nx*ny)
  d=0.001*randomn(seed,nx*ny)
  d=d+dark_value; dark of 3.6 e/pixel
  dark= reform(d,[nx,ny]); in e/pixel
  if n_cosmic ne 0 then ray=cosmic(n_cosmic,60000,nx,ny)
  
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
dark=(t_dark/5)

;creating background frame
backg=cute_background(wavelength=wavelength,token=key, file_bg=file_bg, file_sol=file_sol, file_zod=file_zod)

if(key eq 'default') then token='Default: From STIS backgrounds'
if(key eq'calc') then token= 'Calculated from Zodiacla data'
scaled_spectrum=backg*aeff

;adding background to create an image array
background=dblarr(nx,ny)

for k=0,ny-1 do begin
  background(*,k)=scaled_spectrum
endfor
background=background*exptime*G
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

;generating flats
for k=0,4 do begin
  ;flat
  ;generate flatframe
  ;creating aflat frame
  bias=r_noise*randomn(seed,nx*ny)+b_value
  bias= reform(bias,[nx,ny])
  d=dblarr(nx*ny)
  d=randomn(seed,nx*ny)
  d=d+dark_value; dark of 3.6 e/s/pixel
  dark= reform(d,[nx,ny]); in e/pixel
  flat=5*randomn(seed,nx*ny)+50000L
  flat=reform(flat,[nx,ny])
  if n_cosmic ne 0 then ray=cosmic(n_cosmic,60000,nx,ny)
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
flat=(t_flat/5)

im=make_array(nx,ny, /DOUBLE, value=0.0)
;calculating the scaling value based on transit parameters
if (tag_exist(infile,'transit') eq 0)then begin 
  transit_lightcurve,infile,scale,new_time,mid_time 
endif else begin
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
    
    for m=1, rows-1 do begin
      scale_value=data[*,m]
      scaling=interpol(scale_value,in_wave,wavelength)
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
    endfor
    ;scale=interpol(sacel1,in_wave,ccd_wave) 
  endelse
endelse  
m_t=mid_time/86400.0
Tmid=t+m_t[0]
if (tag_exist(infile,'transit_mid_time') eq 1)then Tmid=double(mid_time)  
;define header of image
;

writecol,file_out+'input_light curve.txt',new_time,scale
t0=t
for k=18, 19 do begin ;n_elements(scale[0,*])-1
;  inc=k*exptime/86400.0
;  t=t0+inc
  if (tag_exist(infile,'transit_mid_time') eq 1)then t=new_time[k] 
  scaled_ccd_count=scale[*,k]*ccd_count
  
  ;defien position on detector
  im=cute_specspread(file_pos,file_spread,scaled_ccd_count,nx,ny,s_pos,0,yshift,pl_scale)
  ;mwrfits,im,file_out+'1sec_image_wtc'+String(k, Format='(I05)') +'.fits',hdr1,/create
  im_wtc=im*exptime*G
  ;adding backgroubd stars
  im_t=im+im_bs
  delta_x=xjitter
  delta_y=yjitter
  delta_z=zjitter
  if infile.jitter_sim eq 1 then im_j=jitter(im_t,delta_x,delta_y,delta_z,exptime) else if infile.jitter_sim eq 0 then im_j=im*exptime

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
 

  
  mwrfits,im_wtc,file_out+'image_wtc'+String(k, Format='(I05)') +'.fits',hdr1,/create
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
  im_n= im_j+sqrt(im_j)
  noi=sqrt(im_j)
    
  im_nb=im_n+background
  ;defining flat
  flat=0.0001*randomn(seed,nx*ny)+1
  flat=reform(flat,[nx,ny])
  ;im_nbf=im_nb/flat ;redundant now
  ;defining bias
  bias=r_noise*randomn(seed,nx*ny)+b_value
  bias= reform(bias,[nx,ny])
  ;im_nbfb=im_nbf+bias
  
  ;defining dark
  d=dblarr(nx*ny)
  d=0.001*randomn(seed,nx*ny)
  d=d+dark_value; dark of 3.6 e/pixel
  dark= reform(d,[nx,ny])
  
  im_nbfbd=(im_nb+(dark*exptime)+bias)*flat
  ;;to counts where G is the gain
  im_final=im_nbfbd*G
  
  ;cosmic ray generation
  ;print,k
  if n_cosmic ne 0 then ray=cosmic(n_cosmic,60000,nx,ny)
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
end