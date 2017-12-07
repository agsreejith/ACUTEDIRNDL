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


pro cutedrndl

;  #Configuration file for the data simulator
;  #
;  #Location of input files
;  input_file      = 'D:\cutedrndl\'
;  #
;  #Stellar Parameters
;  #The parameters for the input star: Temperature in K, Radius in solar radii and V-Magnitude
;  #
;  stellar_temperature = 4100
;  stellar_radius    = 1
;  stellar_magnitude = 10
;  #
;  #Observing Parameters, RA and Dec in degrees, exptime in seconds
;  Ra          = 0
;  Dec         = 0
;  JD          = 0
;  exptime       = 300
;  #
;  #Instrument Parameters :
;  #Position of the star in slit : 3 options avaliable (0',±4',±8' from slit center), input 0,4 or 8
;  #Spectral resolution of the spectrograph in Å.
;  #
;  slit_position   = 0
;  spectral_resolution = 0.8
;  #
;  #Wavelength file location for CUTE
;  wave_file     = 'extra\wavelength.txt'
;  #
;  #Effective area file location for CUTE
;  file_eff='extra\eff_area.txt'
;  #
;  #Spectrum parameter file location for CUTE
;  #Specifies the position on the detector, the slope and y intercept of the the spectrum.
;  #Default file has 3*3 array with each row coresponding to the coresponding position on the slit from center(0)
;  #
;  file_pos='extra\position.txt'
;  #Specify the spread of the spectrum in cross dispresion axis
;  #
;  file_spread='extra\cuteback.txt'
;  #
;  #Background information
;  #Select type of background: Options default or calc,
;  #default: STIS high background values. calc: Calculated from zodiacal light data
;  type_bckg= 'default'
;  #
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
;  #Output location
;  file_out='output\'


;Reading from configuration file

infile=gm_read_textstructure("cutedrndl_parameters.txt");edit the configuration information above and comment this line to run without parameter file.

file=infile.input_file
t_star=infile.stellar_temperature
r_star=infile.stellar_radius
M_star=infile.stellar_magnitude
s_pos=infile.slit_position
fwhm=infile.spectral_resolution
file_wave=infile.wave_file
file_eff=infile.file_eff
file_pos=infile.file_pos
file_spread=infile.file_spread
key=infile.type_bckg
file_bg=infile.bckg_default
file_sol=infile.file_sol
file_zod=infile.file_zod
file_out=infile.file_out
r_noise=infile.readout_noise
exptime=infile.exptime
b_value=infile.bias_value
dark_value=infile.dark_value
G=infile.ccd_gain
ra=infile.Ra
dec=infile.Dec
jd=infile.JD

;get system time
t=systime(/JULIAN, /UTC)
;convert stellar parameters to photons at earth
photons=cute_photons(file,t_star,r_star,M_star)  
photons_star=photons[1,*]
wave=photons[0,*]

photons_star=reform(photons_star,n_elements(photons_star))
wave=reform(wave,n_elements(wave))

;wavelength resolution of CUTE
if s_pos eq 0 then fwhm=fwhm
if s_pos eq 4 then fwhm=fwhm*3/2
if s_pos eq 8 then fwhm=fwhm*2  

hwhm=fwhm/2
;convolution with instrument response
smoothedflux=gaussbroad(wave,photons_star,hwhm)
;smoothedflux = gaussfold(wave, photons_star, fwhm, LAMMIN=2518.50, LAMMAX=3331.10)

w_length=file_lines(file_wave)
wavelength=dblarr(w_length)
ccd_flux=dblarr(2048)
ccd_wave=dblarr(2048)
openr,1,file_wave
readf,1,wavelength
close,1

;interpolate and trim to detector size
ccd_flux=interpol(smoothedflux,wave,wavelength,/SPLINE); linear interpolation
ccd_wave=wavelength

;read effective area and interpolate to cute wavelength
length=file_lines(file_eff)
eff_area=dblarr(11,length)
openr,1,file_eff
readf,1,eff_area
close,1
aeff=interpol(eff_area[4,*],eff_area[0,*],ccd_wave,/SPLINE); linear interpolation

ccd_count=dblarr(2048) 
;Effective area/QE
ccd_count=ccd_flux*aeff

;defien position on detector
im=cute_specspread(file_pos,file_spread,ccd_count,2048,515,s_pos)

im=im*exptime
;background stuff

backg=cute_background(wavelength=wavelength,token=key, file_bg=file_bg, file_sol=file_sol, file_zod=file_zod)

if(key eq 'default') then token='Default: From STIS backgrounds'
if(key eq'calc') then token= 'Calculated from Zodiacla data'
scaled_spectrum=backg*aeff

;adding background to create an image array
background=dblarr(2048,515)

for k=0,514 do begin
  background(*,k)=scaled_spectrum
endfor
background=background*exptime

sxaddpar, hdr1, 'Time in JD', t
sxaddpar, hdr1, 'Star Temperature', t_star
sxaddpar, hdr1, 'Star Radius wrt Rsun', r_star
sxaddpar, hdr1, 'Magnitude ', M_star
sxaddpar, hdr1, 'Slit Position', s_pos
sxaddpar, hdr1, 'Spectral resolution', fwhm
sxaddpar, hdr1, 'Exposure Time in seconds', exptime
;write the fits file
mwrfits,im,file_out+'image_wtc.fits',hdr1,/create

;background header
sxaddpar, bhdr, 'Time in JD', t
sxaddpar, bhdr, 'RA of star (set as 0 is default option is chosen)', ra
sxaddpar, bhdr, 'Dec of star (set as 0 is default option is chosen)', dec
sxaddpar, bhdr, 'Julian date (set as 0 is default option is chosen)', jd
sxaddpar, bhdr, 'Type of background', token
sxaddpar, bhdr, 'Spectral resolution', fwhm
sxaddpar, bhdr, 'Exposure Time in seconds', exptime
;write the fits file
mwrfits,background,file_out+'\background.fits',bhdr,/create

;readoutnoise

;add Photon noise
for i=0,2047 do begin
  for j=0,514 do begin
    im[i,j]= im[i,j]+sqrt(im[i,j])
  endfor
endfor

im=im+background


;flat
;generate flatframe
;creating aflat frame

flat=0.1*randomn(seed,1054720)+1
flat=reform(flat,[2048,515])
sxaddpar, fhdr, 'Time in JD', t
sxaddpar, fhdr, 'Mean', mean(flat)
sxaddpar, fhdr, 'Median ', median(flat)
sxaddpar, fhdr, 'Max value', max(flat)
sxaddpar, fhdr, 'Min Value', min(flat)

;write the fits file
mwrfits,flat,file_out+'flat.fits',fhdr,/create
im=im/flat
;bias
;generate bias

bias=5*randomn(seed,1054720)+r_noise+b_value
bias= reform(bias,[2048,515])
sxaddpar, bihdr, 'Time in JD', t
sxaddpar, bihdr, 'Readout noise', r_noise
sxaddpar, bihdr, 'Bias value', b_value
sxaddpar, bihdr, 'Mean', mean(bias)
sxaddpar, bihdr, 'Median ', median(bias)
sxaddpar, bihdr, 'Max value', max(bias)
sxaddpar, bihdr, 'Min Value', min(bias)


;write the fits file
mwrfits,bias,file_out+'\bias.fits',bihdr,/create

im_b=im+bias

;dark
d=dblarr(1054720)
d=randomn(seed,1054720)
d=d+dark_value; dark of 3.6 e/pixel
dark= reform(d,[2048,515]); in e/pixel

im_f=im_b+dark

;writing dark header
sxaddpar, dhdr, 'Time in JD', t
sxaddpar, dhdr, 'Input dark value', dark_value
sxaddpar, dhdr, 'Mean', mean(dark)
sxaddpar, dhdr, 'Median ', median(dark)
sxaddpar, dhdr, 'Max value', max(dark)
sxaddpar, dhdr, 'Min Value', min(dark)
sxaddpar, dhdr, 'Exposure Time in seconds', exptime
;
;;write the fits file
mwrfits,dark,file_out+'\dark.fits',dhdr,/create

;;to counts where G is the gain

im_f=im_f*G
;
;define header
;
sxaddpar, hdr, 'Time in JD', t
sxaddpar, hdr, 'Star Temperature', t_star
sxaddpar, hdr, 'Star Radius wrt Rsun', r_star
sxaddpar, hdr, 'Magnitude ', M_star
sxaddpar, hdr, 'Slit Position', s_pos
sxaddpar, hdr, 'Spectral resolution', fwhm
sxaddpar, hdr, 'Exposure Time in seconds', exptime
sxaddpar, hdr, 'COMMENT', 'Interstellar Extinction not implimented'
;write the fits file
mwrfits,im_f,file_out+'\raw_image.fits',hdr,/create

print,'Exitng the simulator'

;FITS files creation
end
