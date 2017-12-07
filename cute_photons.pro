; NAME:
;      CUTE_PHOTONS
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

function cute_photons, in_file, t_star, r_star, m_star

;Constants
  R_sun=6.957d10          ; in cm
  L_sun=3.828d33


  t=fix(t_star)
  file=in_file+'\models\t'+String(t, Format='(I05)') +'g4.4\model.flx' ;assuming folder named named by their temperature
  length=file_lines(file)
  fdata=dblarr(3,length)
  openr,1,file
  readf,1,fdata
  close,1
  ;definitions
  r_star=r_star*R_sun
  flux=fdata
  flux[1,*]=3d18*flux[1,*]/(flux[0,*]^2) ;convert to ergs/s/cm2/A
  flux[2,*]=3d18*flux[2,*]/(flux[0,*]^2) ;convert to ergs/s/cm2/A
  flux[1,*]=flux[1,*]*4*!pi*(r_star^2)
  flux[2,*]=flux[2,*]*4*!pi*(r_star^2)
  data=cute_lca(flux,0.3,0.3,10,50)
  wave=data[0,*]
  flux_d=data[1,*]; in erg sec^-1 cm^-2 Hz^-1
  continium=data[2,*]
  photons=dblarr(2,length)
  ;Luminosity and Bolometric magnitude
  L_star=flux_d
  L_bol=total(L_star)
  M_bstar=4.8-2.5*alog10(L_bol/L_sun)
  ;calculate distance to the star


  ;bolometric correction from Flower 1996
  if alog10(t_star) le 3.7 then begin
    BC=-210.13793+0.19596489*t_star-7.4465325e-05*t_star^2+1.4337726e-08*t_star^3- $
      1.3955426e-12*t_star^4+5.4925758e-17*t_star^5
  endif
  if alog10(t_star) ge 4.0 then begin
    BC=4.1940953-0.00070441042*t_star+3.4516521e-08*t_star^2-9.5565244e-13*t_star^3+ $
      1.2790825e-17*t_star^4-6.4741275e-23*t_star^5
  endif
  if alog10(t_star) lt 4.0 and alog10(t_star) gt 3.7 then begin
    BC=-29.325541+0.018052720*t_star-4.4823439e-06*t_star^2.+5.5894085e-10*t_star^3.- $
      3.4753865e-14*t_star^4.+8.5372998e-19*t_star^5
  endif

  M_vstar=M_bstar-BC
  ;calculate E(B-V)
  ;find distance from luminosity and magnitude
  
  ;convert ra dec to
  extn=extinction(glon,glat,distance)
  
  scale=2.512^(M_star-M_vstar) ;to verify
  flux_e=flux_d*scale ;flux at earth

  ;Useful defs
  ;1 Jy = 10^-23 erg sec^-1 cm^-2 Hz^-1
  ;1 Jy = 1.51e7 photons sec^-1 m^-2 (dlambda/lambda)^-1

  print,n_elements(wave)
  photons_star=dblarr(n_elements(wave))

  ;convert to photons
  for i=0,n_elements(wave)-1 do begin
    photons_star[i]=flux_e[i]*5.03e7*wave[i]
  endfor

photons[0,*]=wave
photons[1,*]=photons_star

return,photons
end
