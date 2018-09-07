
pro extinction,glong,glat,distance,av,ebv
;Interstellar Extinction in the Galaxy (Amores & L�pine - 2004)
;This program corresponds to the Axysimetric Model (Model A)
;If you have any difficulty, sugestion or comments, please contact:
;jacques@astro.iag.usp.br     or     amores@astro.iag.usp.br
;You enter longitude, latitude and distance of a point in the Galaxy and get extinction

close,/all
r0=7.5 ;adopted distance of the Galactic center
conv=!pi/180.

step  = 0.05     ;steps of the gas density integration to obtain column density, in pc
;glong=100.0        ;galactic longitude; an arbitrary value given here
;glat=0.0            ;galactic latitude
;dist = 20.0       ;distance of the point to which we will calculate the extinction in kpc

;print,'Interstellar Extinction in the Galaxy (Amores & L�pine - 2005, AJ, 130, 679)'
;
;read,glong,glat,PROMPT='Give the galactic longitude and latitude (Degrees,Degrees)....:  '
;read,dist,PROMPT='Distance [kpc](positive value)...�

dist=distance
nstep=long64(dist/step)

if nstep eq 0 then nstep = 1

;computes  trigonometric functions only once
yproj=cos(glong*conv)
xproj=sin(glong*conv)
bproj=sin(glat*conv)
dproj=cos(glat*conv)
av=0.0                  ;for the integration of the colunar density

;declaring and puting values in the variables. The arrays will contain the
;value of quantities like galactic radius or gas density for each step along the line-of sight
;if you work with other language you should probably define these quantities in a loop
dis= make_array(nstep,/float,value=0.0)
x  = make_array(nstep,/float,value=0.0)
y  = make_array(nstep,/float,value=0.0)
yy = make_array(nstep,/float,value=0.0)
r  = make_array(nstep,/float,value=0.0)
z  = make_array(nstep,/float,value=0.0)
zCO= make_array(nstep,/float,value=0.0)
zH = make_array(nstep,/float,value=0.0)
ah1= make_array(nstep,/float,value=0.0)
aco= make_array(nstep,/float,value=0.0)
zmet=make_array(nstep,/float,value=0.0)
agas=make_array(nstep,/float,value=0.0)
ipas=findGen(nstep)/1 +1  ; generates an array with a sequence of numbers, used as index for
; distance along line-of-sight
nel=n_elements(ipas)

dis=ipas*step - step
x=(dis*xproj)*dproj
y=dis*yproj*dproj
yy=r0-y
r=sqrt(x*x+yy*yy)
z=dis*bproj

zCO=0.036*exp(0.08*r)     ;H2 scale-height
zH = zco*1.8              ;H1 scale-height (Guilbert 1978)
zc = 0.02                 ;shift takes in to account that the sun is not precisely in the galactic plane

ah1=0.7*exp(-r/7.0-((1.9/r)^2))  ;function that calculates the HI density
aco = 58.*exp(-r/1.20-((3.50/r)^2)) + 240.*exp(-(r^2/0.095)) ; H2 density; last term is for galactic center

ah1[0] = 0.0
aco[0] = 0.0

for i=0,nel-1 do begin
  if r[i] le 1.2 then  zmet[i] = 9.6
  if r[i] gt 1.2 and r[i] le 9.0 then zmet[i] = (r0/r[i])^0.5
  if r[i] gt 9.0 then  zmet[i] = (r0/r[i])^0.1
endfor
; this defines the metallicity correction, see section 3 of the paper


gam1=1.0
gam2=2.0

;See the final tuning (section 4.1) correction factor for interval l=120-200

tune=1.
if glong ge 120 and glong le 200 then tune=2.
agas=gam1*(ah1*zmet*exp(-0.5*((z-zc)/zH)^2))+gam2*aco*exp(-0.5*((z-zc)/zCO)^2)

av=total(agas)*step*3.086*.57*tune

; "total" instruction gives the sum of the array elements
; it is equivaletn to integrate along the line-of-sight. The step is in units of kpc=
;3.08 *10^21 cm and the conversion factor gamma= .53 10^-21 mag cm2

rs = 3.05  ;ratio between total to selective extinction
ebv = av/rs

print,'E(B-V)= ',ebv,' mag; Av= ', av, 'mag'

floating_point_underflow = 32
;status = Check_Math()         ; Get status and reset accumulated math error register.
;IF(status AND NOT floating_point_underflow) NE 0 THEN $
;  Message, 'IDL Check_Math() error: ' + StrTrim(status, 2)

return
end
