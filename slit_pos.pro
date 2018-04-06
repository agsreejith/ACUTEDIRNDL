pro slit_pos,ra,dec,ra_s,dec_s,rot_ang,rect_check,xpos,ypos,pos,width,length
close,/all
px=dblarr(4)
py=dblarr(4)
if rot_ang lt 0 then pos_ang=360+rot_ang else pos_ang=rot_ang
theta=pos_ang
theta=theta*!DtoR

w=double(width)/2
l=double(length)/2
l=l/60.0
; calculating both angles

pos=double(pos)
ph1=atan(w/((l-pos)*60))
ph2=atan(w/((l+pos)*60))

R1=(sqrt((w^2)+(((l-pos)*60)^2)))/3600
R2=(sqrt((w^2)+(((l+pos)*60)^2)))/3600

ph = 1.909152433*!DtoR
R = 600.3332408/3600


px[0]=dec+R1*cos(theta+ph1)
px[1]=dec+R1*cos(theta-ph1)
px[2]=dec+R2*cos(!pi+theta+ph2)
px[3]=dec+R2*cos(!pi+theta-ph2)


py[0]=ra+R1*sin(theta+ph1)
py[1]=ra+R1*sin(theta-ph1)
py[2]=ra+R2*sin(!pi+theta+ph2)
py[3]=ra+R2*sin(!pi+theta-ph2)

xpos=0
ypos=0
object = Obj_New('IDLanROI', px, py)
rect_check = object->ContainsPoints(dec_s, ra_s)
Result = object->ComputeGeometry( AREA=area1, CENTROID=centroid1 , PERIMETER=peremeter1)
Obj_Destroy, object
;t2= inside(dec,ra,px,py)

if rect_check eq 1 then begin
  d=sqrt((ra_s-ra)^2+(dec_s-dec)^2)
  alpha = asin((ra_s-ra)/d)-theta
  x=d*sin(alpha)
  y=d*cos(alpha)
  ;print,'d and alpha',ra_s,dec_s,d,alpha
  ypos=y
  xpos=x
;  if rot_ang lt 0 then begin
;    ypos=-y
;    xpos=-x
;  endif
;;  print,ra_s,dec_s
;  print,px
;  print,py
;stop
endif
end