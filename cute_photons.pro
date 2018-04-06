; NAME:
;      CUTE_PHOTONS
;
; PURPOSE:
;      Returns counts in photons/s/cm2 for a star of certain magnitude, radius and temperature ;
; CALLING SEQUENCE:
;      Result = CUTE_PHOTONS(strct)
;
; INPUTS:
;      strct = Structure containing parameter values 
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
;      extinction added 13.12.2017 by A. G. Sreejith

function cute_photons, strct 
close,/all
;Constants
  R_sun=6.957d10          ; in cm
  L_sun=3.828d33

  if (tag_exist(infile,'flux_file') eq 1)then begin
    flux_file=infile.flux_file
    length=file_lines(flux_file)
    fldata=dblarr(2,length)
    openr,1,flux_file
    readf,1,fldata
    close,1
    wave=dblarr(length)
    flux_n=dblarr(length)
    photons=dblarr(length)
  endif else begin
    t_star=strct.stellar_temperature
    t=t_star
    stellar_file=strct.stellar_params
    FMT = 'A,L,F,F'
    readcol,stellar_file,F=FMT,Sp,Teff,BV,Radius
    loc = VALUE_LOCATE(Teff,t)
    if (tag_exist(strct,'stype') eq 0)then sptype = Sp[loc] else sptype = strct.stype
    if (tag_exist(strct,'bv') eq 0) then bv = BV[loc] else bv = strct.bv
    if (tag_exist(strct,'stellar_radius') eq 0) then r_star = Radius[loc] else r_star=strct.stellar_radius
    
    in_file=strct.input_file
    t_star=strct.stellar_temperature
    V_mag=strct.stellar_magnitude
    ra=strct.Ra
    dec=strct.Dec
    jd=strct.JD
    lca_strct={magnitude:V_mag,radius:r_star,temperature:t_star,Ra:ra,Dec:dec,stype:sptype,BV:bv,logr:strct.logr}
    r_star=double(r_star)
    t=fix(t_star)
    
    file=in_file+'\models\t'+String(t, Format='(I05)') +'g4.4\model.flx' ;assuming folder named by their temperature and file are named as model.flx
    if file_test(file) ne 1 then message,'Invalid input in stellar temperature, Please input temperature between 2300K and 10000K in steps of 100'
    length=file_lines(file)
    fdata=dblarr(3,length)
    openr,1,file
    readf,1,fdata
    close,1
    ;definitions
    flux1=dblarr(3,length)
    r_star=double(r_star*R_sun)
    flux=fdata
    flux1[1,*]=(3d18*fdata[1,*])/(fdata[0,*]*fdata[0,*]) ;convert to ergs/cm2/s/A
    flux1[2,*]=(3d18*fdata[2,*])/(fdata[0,*]*fdata[0,*]) ;convert to ergs/cm2/s/A
    flux[1,*]=flux1[1,*]*4*!pi*(r_star^2)  ;convert to ergs/s/A
    flux[2,*]=flux1[2,*]*4*!pi*(r_star^2)  ;convert to ergs/s/A
    t=double(t_star)
    t4=0.0D
    t4=t^4
    r2=r_star^2
    stepahs=5.6704d-5
    L = stepahs*4*!pi*r2*t4
 
    data=cute_lca(flux,lca_strct,0.257,0.288);line core emission 
    wave=data[0,*]
    flux_d=data[1,*]; in erg/s/A
    continium=data[2,*]
    photons=dblarr(2,length)
    ;Luminosity and Bolometric magnitude
    L_star=flux_d
    L_bol=total(L_star)
    M_bstar=4.8-2.5*alog10(L_bol/L_sun)
    ;back to flux density
    ;flux_d=flux_d/(4*!pi*(r_star^2))
    ;writecol,'flux_'+String(t, Format='(I05)') +'.txt',wave,flux[1,*],flux_d

    L2= fdata[1,*]*!pi*(r_star^2)
    L2_bol=total(L2)
    T=double(t_star)
    ;bolometric correction from Flower 1996
    if alog10(T) le 3.7 then begin
      BC=-210.13793+0.19596489*T-7.4465325e-05*T^2+1.4337726e-08*T^3- $
      1.3955426e-12*T^4+5.4925758e-17*T^5
    endif
    if alog10(T) ge 4.0 then begin
      BC=4.1940953-0.00070441042*T+3.4516521e-08*T^2-9.5565244e-13*T^3+ $
      1.2790825e-17*T^4-6.4741275e-23*T^5
    endif
    if alog10(T) lt 4.0 and alog10(T) gt 3.7 then begin
      BC=-29.325541+0.018052720*T-4.4823439e-06*T^2.+5.5894085e-10*T^3.- $
      3.4753865e-14*T^4.+8.5372998e-19*T^5
    endif

    M_vstar=M_bstar-BC
    
    ;find distance from luminosity and magnitude
    d=10*10^(0.2*(V_mag-M_vstar))
    ;scale=2.512^(M_star-M_vstar) ;to verify
    ;convert ra dec to galactic  
    euler,dec,ra,glon,glat,1
  
    ;  da=mrdfits('extra\vmag.fits',0,h)
    ;  f_wave=da[0,*]
    ;  f_flx=da[1,*]
    ;  fl_wave=reform(f_wave,n_elements(f_wave))
    ;  fl_flx=reform(f_flx,n_elements(f_flx))
    ;  wavelength=reform(wave,n_elements(wave))
    ;  wmin = min(fl_wave, max = wmax)
    ;  g = where((wavelength GE wmin) and (wavelength LE wmax), Ng)
    ;  if Ng EQ 0 then begin
    ;    message,/CON, $
    ;      'ERROR - Supplied wavelengths outside of filter response curve'
    ;    return,-1
    ;  endif
    ;
    ;  wgrid = [fl_wave,wavelength(g)]
    ;
    ;  wgrid = wgrid(sort(wgrid))
    ;  wgrid = wgrid(uniq(wgrid))
    ;flux_el=reform(flux_e,n_elements(flux_e))
    ;L_starl=reform(L_star,n_elements(L_star))
    ;  linterp,wavelength,flux_el,wgrid,fluxg
    ;  linterp,fl_wave,fl_flx,wgrid,frelg
    ;  Lv = fluxg*frelg              ;Effective flux through filter
    ;  Ls=tsum(wgrid,Lv)
    ;  linterp,wavelength,L_starl,wgrid,fluxg2
    ;  linterp,fl_wave,fl_flx,wgrid,frelg2
    ;  Bv = fluxg2*frelg2
    ;  Bs=tsum(wgrid,Bv)
    ;
    ;
    ;d=sqrt(Bs/(4*!pi*Ls))
    
    ;calculate E(B-V)
    extinction,glon,glat,d,av,ebv ;get extinction based on coordinates and distance 
    
    if (tag_exist(strct,'mg_col') eq 1)then nmg2 = strct.mg_col else begin
       nh=5.8d21*ebv
       ;The Mg2 column density is
       fractionMg2 = 0.825 ;   (Frisch & Slavin 2003; this is the fraction of Mg in the ISM that is singly ionised)
       Mg_abn      = -5.33 ;   (Frisch & Slavin 2003; this is the ISM abundance of Mg)
       nmg2        = alog10(nh*fractionMg2*10.^Mg_abn)
    endelse
   
    
    wave_data=dblarr(length)
    flux_f=dblarr(3,length)
    flux_f[0,*]=wave
    flux_f[1,*]=flux_d
    flux_f[2,*]=continium
    flux_data=cute_ism_abs(data,nmg2)
    flux_di=flux_data[1,*]
    
    flux_e=flux_di/(4*!pi*(d*(3.086e+18)^2)) ;flux at earth
    escale=1
    ;escale=10^((av)/2.5)
    ebv=-1*ebv
    fm_unred,wave,flux_e,ebv,flux_n ;flux at earth
    ;flux_n=flux_e/escale

    ;Useful defs
    ;1 Jy = 10^-23 erg sec^-1 cm^-2 Hz^-1
    ;1 Jy = 1.51e7 photons sec^-1 m^-2 (dlambda/lambda)^-1
  endelse
    photons_star=dblarr(n_elements(wave))
    ;convert to photons
    photons_star=flux_n*5.03e7*wave ;from ergs/s/cm2/A to photons/s/cm2/A
    photons[0,*]=wave
    photons[1,*]=photons_star
 
    return,photons
  
end