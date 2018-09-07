pro cute_field,ra_t,dec_t,s_ra,s_dec,rot_ang,file_out,pos,v_mag,M_star,width,length
close,/all
!x.style=3        ; Set explicit limits
!y.style=3
px=dblarr(5)
py=dblarr(5)
ra=double(ra_t)
dec=double(dec_t)
if rot_ang lt 0 then theta=360+rot_ang else theta=rot_ang
theta=theta*!DtoR
; calculating both angles
w=double(width)/2
l=double(length)/2
l=l/60.0
pos=double(pos)
ph1=atan(w/((l-pos)*60))
ph2=atan(w/((l+pos)*60))

R1=(sqrt((w^2)+(((l-pos)*60)^2)))/3600
R2=(sqrt((w^2)+(((l+pos)*60)^2)))/3600


;ph = 1.909152433*!DtoR
;R = 600.3332408/3600

px[0]=dec+R1*cos(theta+ph1)
px[1]=dec+R1*cos(theta-ph1)
px[2]=dec+R2*cos(!pi+theta+ph2)
px[3]=dec+R2*cos(!pi+theta-ph2)
px[4]=dec+R1*cos(theta+ph1)

py[0]=ra+R1*sin(theta+ph1)
py[1]=ra+R1*sin(theta-ph1)
py[2]=ra+R2*sin(!pi+theta+ph2)
py[3]=ra+R2*sin(!pi+theta-ph2)
py[4]=ra+R1*sin(theta+ph1)

t_d=make_array(3,/DOUBLE, value=dec)
t_r=make_array(3,/DOUBLE, value=ra)


req = where(v_mag lt 20)
req2= where(v_mag eq 22)
n_v=v_mag(req)
n_ra=s_ra(req)
n_dec=s_dec(req)
nq_v=v_mag(req2)
nq_ra=s_ra(req2)
nq_dec=s_dec(req2)
;;intensity color
;colors = cgScaleVector(Findgen(N_Elements(n_v)), Min(n_v), Max(n_v))
;clr = Value_Locate(colors, n_v)
;clr = Byte(Round(cgScaleVector(clr, 0, 255)))
;cgLoadCT, 33
;cgWindow
cgDisplay, 1400,1400
;!p.position=[0.15,0.15,0.99,0.99]
m=max(v_mag)
mi=min(v_mag)
range=m-mi
scale=range/100

value=20/n_v
!x.style=3        ; Set explicit limits
!y.style=3
;the scatterplot
theta = 2*!pi * findgen(21)/20
usersym, 1*cos(theta), 1*sin(theta) $
  , fill=0, thick=0
  
cgplot, n_ra, n_dec,background=cgColor('white'),yrange=[dec-0.20000001,dec+0.1999999],xrange=[ra+0.199999,ra-0.2000001],color=cgcolor('black'),xtitle='right ascension [degree]', $
  ytitle='declination [degree]',title='Ra:'+String(ra, format='(F9.5)')+', Dec:'+String(dec, format='(F9.5)')+', Pos angle:'+String(rot_ang, format='(F5.1)')+', Slit pos:'+String(pos, format='(F3.1)'),$
  psym=3,symsize=1,charsize=3,charthick=2.5
  size_val=20/double(M_star)
  cgoplot,ra,dec,psym=16,color=cgColor('medium grey'),symsize=size_val
for i=0, n_elements(n_dec)-1 do begin
  p_d=make_array(3,/DOUBLE, value=n_dec[i])
  p_r=make_array(3,/DOUBLE, value=n_ra[i])
  cgoplot,n_ra[i],n_dec[i],color=cgcolor('black'), psym=9,symsize=value[i]
endfor

cgoplot,nq_ra,nq_dec,psym=7,symsize=0.5 
theta = 2*!pi * findgen(21)/20
usersym, 1.0*cos(theta), 1.0*sin(theta) $
  , fill=1, thick=1
;cgplot,s_dec,s_ra,xrange=[dec-0.2,dec+0.2],yrange=[ra-0.2,ra+0.2],psym=16,ytitle='RA',xtitle='Declination',background=cgColor('white'),color=cgColor('black'),thick=1.5,xthick=1.2,ythick=1.2,charsize=2,charthick=1.5
;cgoplot,ra,dec,psym=16,color=cgColor('medium grey'),symsize=size_val
;cgColorbar, Divisions=10, Format='(I3)', Range=[min(n_v),max(n_v)], /AddCMD , Title='V_mag'

cgoplot, n_ra, n_dec,psym=3,symsize=1
cgPolygon,py,px,color=cgColor('blue')
cgArrow, ra-0.19, dec-0.19, ra-0.14,dec-0.19, /Solid,/DATA
cgArrow, ra-0.19, dec-0.19, ra-0.19, dec-0.14, /Solid,/DATA
cgText,ra-0.19, dec-0.135,'N',CHARSIZE=1.5,ALIGNMENT=0.5,/DATA
cgText,ra-0.135, dec-0.19,'E',CHARSIZE=1.5,ALIGNMENT=0.5,/DATA
;xyouts,100,50,'N',CHARSIZE=10
dia=string(M_star)
;dia2=string(20)
cgLegend, Title=['V:20'], PSym=[9], SymSize=[1], Color=['black'], Location=[0.84, 0.88], Length=0.0, VSpace=1.0, /Box,/Background, BG_Color='white'
cgLegend, Title=['V:'+dia], PSym=[9], SymSize=[size_val], Color=['black'], Location=[0.14, 0.88], Length=0.0, VSpace=1.0, /Box,/Background, BG_Color='white'

file_var=fix(rot_ang)
write_png,file_out+'cute_field_'+String(file_var, Format='(I03)')+'.png',TVRD(/TRUE)
stop
end