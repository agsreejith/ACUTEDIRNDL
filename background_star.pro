; NAME:
;      BACKGROUND_STAR
;
; PURPOSE:
;      Creates image array for background stars in CUTE spectrograph;
; CALLING SEQUENCE:
;      Result = BACKGROUND_STAR(structure)
;
; INPUTS:
;      structure: Structure containing parameter values
;
; OUTPUT:
;      Returns image array incuding background stars.
; REQUIRES:
;      The parameters of the star is required.
;
; PROCEDURE:
;      Create background star image for CUTE spectrograph for a particular observation;


function background_star,infile
close,/all
t2=systime(/SECONDS)

;old method of aladin retrivel
;aladin_link='http://simbad.u-strasbg.fr/simbad/sim-coo?output.format=ASCII&list.spsel=on&list.mtsel=on&obj.spsel=on&obj.fluxsel=on&v=on&u=off&otypedisp=S&Coord=12%2030%20%2b10%2020&Radius=20&Radius.unit=arcmin'
;
;a=webget(aladin_link)
;openw,1,'aladin.txt'
;printf,1,a.Text
;close,1
;aladin_template={ VERSION: 1.0000000,$
;                  DATASTART: 13,$
;                  DELIMITER: '|',$
;                  MISSINGVALUE: !VALUES.F_NAN,$
;                  COMMENTSYMBOL: "=",$
;                  FIELDCOUNT: 14,$
;                  FIELDTYPES: [3, 4, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 3, 3],$
;                  FIELDNAMES: [ "num", "distance", "identifier", "type", "coordinate", "u", "b", "v", "r", "i", "spec", "morph", "bib", "FIELD14"],$
;                  FIELDLOCATIONS: [0, 4, 15, 48, 62, 90, 97, 104, 111, 118, 125, 141, 173, 178],$
;                  FIELDGROUPS: [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13]$
;}
;data = READ_ASCII('aladin.txt', TEMPLATE=aladin_template)
;print,data.coordinate[2]
;print,n_elements(data.coordinate)
;for i=0, n_elements(data.coordinate)-1 do begin
;  t=strsplit(data.coordinate[i], /EXTRACT)
;  print,i
;  ra[i]=double(t[0])+(double(t[1])/60)+(double(t[2])/3600)
;  dec[i]=double(t[3])+(double(t[4])/60)+(double(t[5])/3600)
;  ;print,ra,dec
;endfor

;ra=strct.Ra
;dec=strct.Dec
;jd=strct.JD
;in_file=strct.input_file
;t_star=strct.stellar_temperature
;r_star=strct.stellar_radius
;V_mag=strct.stellar_magnitude


;definitions from parameter file
in_file=infile.input_file
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
ra=infile.Ra
dec=infile.Dec
jd=infile.JD
exptime=infile.exptime
s_wid=infile.slit_width
s_len=infile.slit_length
;jitter_val=infile.jitter
nx=long(infile.x_pixels)
ny=long(infile.y_pixels)
star_position=fix(infile.slit_position)
theta=double(infile.rotation)
pl_scale=double(infile.plate_scale)

diff_value=dblarr(3)
im_str=dblarr(nx,ny)
w_length=file_lines(file_wave)
wavelength=dblarr(w_length)
ccd_flux=dblarr(nx)
ccd_wave=dblarr(nx)
openr,1,file_wave
readf,1,wavelength
close,1

if (tag_exist(infile,'bckg_star') eq 0)then begin
  ;get catalogue from VIZIER
  data=QUERYVIZIER('I/305/out',[ra,dec],12,/ALLCOLUMNS)
  tags = Get_Tags(data)
  ;._R .RAJ2000 .DEJ2000 .SRCIDGAIA .RAGAIA .E_RAGAIA .DEGAIA .E_DEGAIA .ORG .NU .EPUCAC .PMRA .E_PMRA .PMDE .E_PMDE .GMAG .F_MAG .RMAG .JMAG .HMAG .KMAG
;print,tags

  ;retriving data from catalog
  s_ra=data.RAJ2000
  s_dec=data.DEJ2000
  vmag=data.vmag
  bmag=data.bmag 
  nmag=data.nmag

  endif else begin
    bckg_str_file=infile.bckg_str
    length=file_lines(bckg_str_file)
    bs_data=dblarr(5,length)
    openr,1,bckg_str_file
    readf,1,bs_data
    close,1
    s_ra=dblarr(length)
    s_dec=dblarr(length)
    sr_temp=dblarr(length)
    sr_mag=dblarr(length)
    sr_bv=dblarr(length)
    s_ra=bs_data[0,*]
    s_dec=bs_data[1,*]
    sr_temp=bs_data[2,*]
    sr_mag=bs_data[3,*]
    sr_bv=bs_data[4,*]
  endelse
;constants for converting from 2MASS to uvbri based on bilir et al 2007   arXiv:0711.4356
;a1 = 1.210 
;a2 = 1.816 
;a3 = 1.896 
;b1 = 1.295
;b2 = 1.035 
;b3 = 1.131
;c1 = -0.046
;c2 = 0.016 
;c3 = -0.004 
;A1 = a2-a1
;A2 = a2-a3
;B1 = b2-b1
;B2 = b2-b3
;C1 = c2-c1
;C2 = c2-c3


;reading the spectral type file
;spectra_file='extra\spectral_type.txt'
;l=file_lines(spectra_file)
;spectral_type=strarr(l)
;temperature=make_array(l,/L64)
;mag=fltarr(l)
;lum=dblarr(l)
;FMT='A,LL,F,D'
;READCOL,spectra_file,F=FMT,spectral_type,temperature,mag,lum
;calculate the counts for each background star 
FMT = 'A,L,F,F'
readcol,infile.stellar_params,F=FMT,Sp,Teff,BV,R


for j=0, n_elements(s_ra)-1 do begin
  slit_pos,ra,dec,s_ra[j],s_dec[j],theta,rect_check,star_pos_x,star_pos_y,s_pos,s_wid,s_len ;checking if star is within the slit and getting the x and y dispesion from the target star
  seperation=sqrt((s_ra[j]-ra)^2+(s_dec[j]-dec)^2)
  positon=seperation*3600
  if positon gt 1 then begin ; checking to remove target star
    if (tag_exist(infile,'bckg_star') eq 0)then begin
      if ((finite(vmag[j],/NaN) eq 0) and (finite(bmag[j],/NaN) eq 0)) then begin  
        v_mag=vmag[j]
        s_bv=bmag[j]-vmag[j]
        loc = VALUE_LOCATE(BV,s_bv)
        s_temp = Teff[loc]
        s_rad = R[loc]  ;stellar radi wrt sun
      endif else begin
        if ((finite(bmag[j],/NaN) eq 0) or (finite(nmag[j],/NaN)) eq 0)then begin
         ;print,bmag[j],nmag[j] 
         if ((bmag[j] le 18.0) or (nmag[j] le 18.0)) then begin 
            vmag[j] = 22.0
          endif else vmag[j]=25.0
        endif else vmag[j]=25
        v_mag=25
        goto, skip
      endelse 
    endif else begin
      v_mag=sr_mag[j]
      s_temp=sr_temp[j]
      loc=VALUE_LOCATE(Teff,s_temp)
      s_rad = R[loc]
    endelse
      
        if rect_check eq 1 then begin 
         ;calculating stellar parameters
        ;add stuff for Line core emission values
        ;print,s_ra[j],s_dec[j],v_mag,star_pos_y*3600,s_rad,s_temp
        
          in_structure={input_file:in_file,stellar_magnitude:v_mag,stellar_radius:s_rad,stellar_temperature:s_temp,Ra:s_ra[j], $
          Dec:s_dec[j],JD:jd,stype:Sp[loc],stellar_params:infile.stellar_params,logr:-4.9}
          if  (tag_exist(infile,'mg_col') eq 1)then struct_add_field,in_structures, mg_col,infile.mg_col  
            wave=wavelength
            bs_photons=cute_photons(in_structure)
            bs_photons_star=bs_photons[1,*]
            bs_wave=bs_photons[0,*]
  
            bs_photons_star=reform(bs_photons_star,n_elements(bs_photons_star))
            bs_wave=reform(bs_wave,n_elements(bs_wave))
            ;print,s_ra[j],s_dec[j],v_mag,star_pos_y*3600
      
            ;assuming the slit is oriented parallel to RA axis with object at the center
            if ((star_pos_y*60) lt 4) then pos = 0   
            if (((star_pos_y*60) lt 8) and ((star_pos_y*60) ge 4))then pos = 4
            if ((star_pos_y*60) ge 8) then pos = 8
                        
            ;wavelength resolution of CUTE
            if (s_pos ge -2 && s_pos le 2) then fwhm=fwhm
            if (s_pos gt 2 && s_pos le 8) then fwhm=fwhm*3/2
            if (s_pos gt 8) then fwhm=fwhm*2
            if (s_pos lt -2 && s_pos ge -8) then fwhm=fwhm*3/2
            if s_pos lt -8 then fwhm=fwhm*2 
            st=value_locate(bs_wave, 2000)
            en=value_locate(bs_wave, 3500)
            bs_wave_new=bs_wave[st:en]
            bs_photons_star_new=bs_photons_star[st:en]
            hwhm=fwhm/2
            ;convolution with instrument response
            bs_smoothedflux=gaussbroad(bs_wave_new,bs_photons_star_new,hwhm)
            x_shift = fix(star_pos_x*(3600/pl_scale)) ;in pixel 
            y_shift = fix((star_pos_y+(star_position/60))*(3600/pl_scale)) 
            new_wave=dblarr(nx+ ABS(x_shift))
            if x_shift lt 0 then begin
              xshift=ABS(x_shift)
              for i=0,ABS(x_shift)-1 do begin
                diff_value[0]=wave[1]-wave[0]
                diff_value[1]=wave[2]-wave[1]
                diff_value[2]=wave[3]-wave[2]
                diff = mean(diff_value)
                new_wave[i]= -diff*(i+1)+wave[0]
              endfor  
              new_wave(xshift:nx-1+xshift)=wave
              bs_flux=interpol(bs_smoothedflux,bs_wave_new,new_wave,/SPLINE);
              for i=0,nx-1 do begin
                ccd_flux[i]=bs_flux[i]
                ccd_wave[i]=new_wave[i]
              endfor
            endif
            if x_shift gt 0 then begin
              xshift=ABS(x_shift)
              new_wave[nx-1]=wave[nx-1]
              for i=nx,nx+x_shift-1 do begin
                diff_value[0]=wave[nx-1]-wave[nx-2]
                diff_value[1]=wave[nx-2]-wave[nx-3]
                diff_value[2]=wave[nx-3]-wave[nx-4]
                diff = mean(diff_value)
                new_wave[i]= diff+new_wave[i-1]
              endfor
              new_wave[0:nx-1]=wave
              bs_flux=interpol(bs_smoothedflux,bs_wave_new,new_wave,/SPLINE);
              for i=0,nx-1 do begin
                ccd_flux[i]=bs_flux[i+x_shift]
                ccd_wave[i]=new_wave[i+x_shift]
              endfor
            endif 
            if x_shift eq 0 then begin
              bs_flux=interpol(bs_smoothedflux,bs_wave_new,wave,/SPLINE);
              ccd_flux=bs_flux
              ccd_wave=wave
            endif
            length=file_lines(file_eff)
            eff_area=dblarr(11,length)
            openr,1,file_eff
            readf,1,eff_area
            close,1
            aeff=interpol(eff_area[4,*],eff_area[0,*],ccd_wave,/SPLINE); linear interpolation
            ccd_count=dblarr(nx)
            ;Effective area/QE
            ccd_count=ccd_flux*aeff
            ;if subs eq 2045 then stop
            ;defien position on detector
            im=cute_specspread(file_pos,file_spread,ccd_count,nx,ny,s_pos,0,y_shift,pl_scale)
            im_str+=im
        endif 
    skip: 
  endif
endfor
writecol,file_out+'background_star.txt',s_ra,s_dec,vmag
cute_field,ra,dec,s_ra,s_dec,theta,file_out,s_pos,vmag,M_star,s_wid,s_len
;jitter is implimented in main code
;delta_x=jitter_val
;delta_y=jitter_val
;delta_z=jitter_val
;t3=systime(/SECONDS)
;print,'jitter2',t3-t2
;if infile.jitter_sim eq 1 then im_jstr=jitter(im_str,delta_x,delta_y,delta_z,exptime) else im_jstr= im_str*exptime 
;t4=systime(/SECONDS)
;print,'jitter2',t4-t3

;return the background for 1 second exposure

return,im_str
end