; NAME:
;      LIMB_DARKENING
;
; PURPOSE:
;      Returns the linear and quadratic limb-darkening coefficient for a star of particular temperature in CUTE bandwidth;
; CALLING SEQUENCE:
;       LIMB_DARKENING,strct,u1,u2
;
; INPUTS:
;      strct = Structure containing parameter values
;
; OUTPUT:
;      u1: linear limb-darkening coefficient for each wavelength
;      u2: quadratic limb-darkening coefficient for each wavelength
; REQUIRES:
;      The input stellar spectrum files stored at exact locations. Uses mpfitfun for fitting
;
; PROCEDURE:
;      Calculate the linear and quadratic limb-darkening coefficient for a star of particular temperature in CUTE bandwidth.
; MODIFICATION HISTORY:
;      created 21.01.2018 by A. G. Sreejith
;      modified 05.02.2018 
function ld_quad, x, p
  return, 1. - p[0]*(1.-x) - p[1]*(1.-x)^2
end
function ld_nl, x, p
  return, 1. - p[0]*(1.-x^(0.5)) - p[1]*(1.-x) - p[2]*(1.-x^(1.5)) - p[3]*(1.-x^2)
end

function ma_ang, x, p
  ;x=r/rs
  u1=p[0]
  u2=p[1]
  return,[1.-u1*(1-sqrt(1-(sin(acos(x)))^2))-u2*(1-sqrt(1-(sin(acos(x)))^2))^2]/(1-u1/3-u2/6)/!pi
end

pro limb_darkening,infile,n_u1,n_u2  
close,/all

;infile=gm_read_textstructure("cutedrndl_parameters.txt")
file_wave=infile.wave_file
file_eff=infile.file_eff
t_star=infile.stellar_temperature
l_file=infile.limb_file
ld_fit=1   ;limb darkening coefficients: 1=quadratic law & 2=non-linear law & 3 Mandel & Agol (2002) coeffs

w_length=file_lines(file_wave)
wavelength=dblarr(w_length)
ccd_flux=dblarr(w_length)
ccd_wave=dblarr(w_length)

openr,1,file_wave
readf,1,wavelength
close,1
ccd_wave=wavelength
;file_eff='extra\eff_area.txt'
length=file_lines(file_eff)
eff_area=dblarr(11,length)
openr,1,file_eff
readf,1,eff_area
close,1
aeff=interpol(eff_area[4,*],eff_area[0,*],ccd_wave,/SPLINE);

aeff=aeff/max(aeff)
bandwidth=ccd_wave[w_length-1]-ccd_wave[0]
bin = bandwidth/16
bin=w_length/16



t=fix(t_star)
phonex_file=l_file+'\lte'+String(t, Format='(I05)')+'-4.50-0.0.PHOENIX-ACES-AGSS-COND-SPECINT-2011.fits'
if file_test( phonex_file) ne 1 then t = t+100
phonex_file=l_file+'\lte'+String(t, Format='(I05)')+'-4.50-0.0.PHOENIX-ACES-AGSS-COND-SPECINT-2011.fits'
mu=mrdfits(phonex_file,1,head)
flux_all=mrdfits(phonex_file,0,head)
;print,head
wl0=sxpar(head,'CRVAL1')
dwl=sxpar(head,'CDELT1')
;print,wl0,dwl
;wl0=hierarch(head,'CRVAL1') & dwl=hierarch(head,'CDELT1')
wl_ld=dblarr(n_elements(flux_all[*,0]))
for i=0,n_elements(wl_ld)-1 do wl_ld[i]=wl0+i*dwl


valid=where(mu ge 0)
u1=dblarr(w_length/bin)
u2=dblarr(w_length/bin)
wl=dblarr(w_length/bin)
jk=0
for j=0L, w_length-1, bin do begin
    filter=dblarr(2,10)
    filter[0,*]=ccd_wave[j:j+9]
    filter[1,*]=aeff[j:j+9]

    flux_new=dblarr(n_elements(filter[0,*]),n_elements(flux_all[0,*]))
    for i=0,n_elements(flux_all[0,*])-1 do begin
      flux_new[*,i]=interpol(flux_all[*,i],wl_ld,filter[0,*])*filter[1,*]
    endfor

   ; mu=(mu-mu[22])/(1.-mu[22])

    flux=dblarr(n_elements(mu))
    for i=0,n_elements(mu)-1 do flux[i]=total(flux_new[*,i])
;    plot,mu,flux/flux[n_elements(mu)-1],thick=2;,/ylog
;    oplot,mu,flux/flux[n_elements(mu)-1],thick=2,psym=2,symsize=2
;vline,0.

    err=mu & for i=0,n_elements(err)-1 do err[i]=1d-15

    ii=where(mu gt 0.07 and flux gt 0.05)

    if ld_fit eq 1 then begin
      start=[1.5,-0.5]
      ld=mpfitfun('ld_quad', mu[ii], flux[ii]/flux[n_elements(mu)-1], err[ii], start) 
      ;oplot,mu,ld_quad(mu,ld),psym=6,thick=2
      
      ;oplot,mu,1. - 0.4952*(1.-mu) - 0.3003*(1.-mu)^2,color=4,thick=2
     
    endif

    if ld_fit eq 3 then begin
      start=[1.0,1.5]
      ld=mpfitfun('ma_ang', mu[ii], flux[ii]/flux[n_elements(mu)-1], err[ii], start)
      
      ;oplot,mu,1. - start[0]*(1.-mu^(0.5)) - start[1]*(1.-mu) - start[2]*(1.-mu^(1.5)) - start[3]*(1.-mu^2),color=4,thick=2
      ;oplot,mu,1. - 0.25*(1.-mu^(0.5)) - 1.02*(1.-mu) + 0.5*(1.-mu^(1.5)) - 0.08*(1.-mu^2),color=4,thick=2
      
      ;print,ld
    endif

    if ld_fit eq 2 then begin
      start=[0.0,1.5,0.0,-0.5]
      ld=mpfitfun('ld_nl', mu[ii], flux[ii]/flux[n_elements(mu)-1], err[ii], start)
      
      ;oplot,mu,1. - start[0]*(1.-mu^(0.5)) - start[1]*(1.-mu) - start[2]*(1.-mu^(1.5)) - start[3]*(1.-mu^2),color=4,thick=2
      ;oplot,mu,1. - 0.25*(1.-mu^(0.5)) - 1.02*(1.-mu) + 0.5*(1.-mu^(1.5)) - 0.08*(1.-mu^2),color=4,thick=2
      
      ;print,ld
    endif
    
    u1[jk]=ld[0]
    u2[jk]=ld[1]
    wavecent=(j+j+9)/2
    wl[jk]=ccd_wave[wavecent]
    ;LIMB DARKENING
    ;ii=where( (wl_ld ge 5500. and wl_ld le 5549.) or (wl_ld ge 5555. and wl_ld le 5602.) or (wl_ld ge 5608. and wl_ld le 6800.) )
    ;flux=dblarr(n_elements(mu))
    ;for i=0,n_elements(mu)-1 do flux[i]=total(flux_all[ii,i])
    ;plot,1.-mu,flux/flux[n_elements(mu)-1]
    ;
    ;err=mu & for i=0,n_elements(err)-1 do err[i]=1d-15
    ;start=[0.3,0.3]
    ;ld_white=mpfitfun('myld', mu[50:77], flux[50:77]/flux[n_elements(mu)-1], err[50:77], start)
    ;oplot,1.-mu,myld(1.-mu,ld_white),color=2
jk++
endfor

n_u1=interpol(u1,wl,ccd_wave)
n_u2=interpol(u2,wl,ccd_wave)

end







