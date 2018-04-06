; NAME:
;      TRANSIT_LIGHTCURVE
;
; PURPOSE:
;      Returns the scaling value for flux at each wavelength during transit in CUTE bandwidth;
; CALLING SEQUENCE:
;       TRANSIT_LIGHTCURVE,strct,scale
;
; INPUTS:
;      strct = Structure containing parameter values
;
; OUTPUT:
;      scale: averaged scaling value (for exposure time) at each wavelength for transit
;      
; REQUIRES:
;      The input stellar spectrum files stored at exact locations, uses exofast_occultquad.pro for transit calculations
;
; PROCEDURE:
;      Calculate the averaged scaling value (for exposure time) at each wavelength for transit;
; MODIFICATION HISTORY:
;      created 22.01.2018 by A. G. Sreejith
;      modified 05.02.2018 

pro transit_lightcurve,infile,scale,new_time,Tmid
close,/all
;infile=gm_read_textstructure("cutedrndl_parameters.txt")
R_sun=6.957d10 
nx      = long(infile.x_pixels)
b       = double(infile.impact_parameter)
Porb    = double(infile.orbital_period)
incli   = double(infile.orbital_inclination)
a       = double(infile.rlt_semi_major_axis)
exptime = fix(infile.exptime)
R_star  = double(infile.stellar_radius)
r_time  = fix(infile.read_time)
time    = dblarr(40001);defining time array with arbitary large size
scaling = dblarr(nx,40001) ; defining a large scaling array
sini    = sin(incli*!pi/180d0)
x       = (dindgen(40001)/40000d0 - 0.5d0)*4d0
z       = sqrt(x^2 + b^2)
a       = a*R_sun*R_star

;b=(a*cos(incli*!pi/180d0))/(R_star*R_sun)

 if (tag_exist(infile,'planet_radius') eq 0)then begin
  file_wave=infile.wave_file
  planet_file=infile.planet_radius_file
  len=file_lines(planet_file)
  data=dblarr(2,len)
  openr,1,planet_file
  readf,1,data
  close,1
  w_length=file_lines(file_wave)
  wavelength=dblarr(w_length)
  openr,1,file_wave
  readf,1,wavelength
  close,1
  pl=interpol(data[1,*],data[0,*],wavelength)
 endif else begin
  p = double(infile.planet_radius)
    pl=make_array(nx,/DOUBLE, value=0.25);modification for  hd209458
 endelse
 ;R_pl    = mean(pl)*R_sun*R_star 
  R_pl=p*R_sun*R_star
 t_dur  = double((Porb/!pi)*asin(sqrt(((R_star*R_sun+R_pl)^2)-(b^2))/a*sini))

  ;for i=0L,n do time[n-i]=Tmid-i*t_step
  ;for i=1L,n do time[n+i]=Tmid+i*t_step
  limb_darkening,infile,u1,u2
  ii=where(z eq b)
 

  for k=0,nx-1 do begin
    exofast_occultquad, z, u1[k], u2[k], pl[k], muo1, muo, d=d
    scaling[k,*]=muo1
  endfor  
    vec1=muo1[0:ii]
    vec2=muo1[ii+1:n_elements(z)-1]
    t1_ar=where(vec1 lt 1.0)
    t2_ar=where(vec2 lt 1.0)
    t1=t1_ar[0]-1
    t2=t2_ar[n_elements(t2_ar)-1]+1
    interval=t2+ii-t1+1
    step=t_dur/interval
    for i=0L, 40000 do time[i]=i*step*86400 ;time in seconds
    writecol,'output\light_curve1.txt',time,muo1
    obs_time=time[n_elements(time)-1]-time[0]
    rounding = obs_time mod (exptime+r_time)
    reminder=exptime+r_time-rounding
    end_time = obs_time-rounding
    steps=value_locate(time, exptime+r_time)
    en= value_locate(time,end_time)
    j=0
    exp_step=value_locate(time,exptime)
    new_time=dblarr((en/steps)-1)
    scale=dblarr(nx,(en/steps)-1)
    for k=0,nx-1 do begin
      jj=0
      for j=0, n_elements(scale[0,*])-1 do begin
        ;print,jj,j,jj+exp_step
        new_impct=scaling[k,jj:jj+exp_step]
        scale[k,j]=mean(new_impct)
        jj=jj+steps
        new_time[j]=time[jj]
        loc=where(x eq 0)
        Tmid=time[loc]
       endfor
    endfor 
if (tag_exist(infile,'transit_mid_time') eq 1)then begin
  
  loc=where(x eq 0)
  Tmid=infile.transit_mid_time
  for i=0L,20000 do time[20000-i]=Tmid-i*step
  for i=1L,20000 do time[20000+i]=Tmid+i*step
  j=0
  for i=0L,(n_elements(scale[0,*])-1)*steps,steps do begin
    new_time[j]=time[i]
    j++
  endfor
endif

end