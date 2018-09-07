pro cute_systematics,infile,t1,t0,im_in

;  f(x) = p1*x^3 + p2*x^2 + p3*x + p4
;  Coefficients (with 95% confidence bounds):
;  p1 =  -3.453e-12  (-6.928e-12, 2.124e-14)
;  p2 =   4.092e-08  (1.734e-10, 8.167e-08)
;  p3 =  -0.0001206  (-0.0002752, 3.391e-05)
;  p4 =         1.1  (0.9112, 1.289)

s=size(im_in)
nx=s[1]
ny=s[2]
im_out=dblarr(nx,ny)
  p1=infile.sym_p1
  p2=infile.sym_p2
  p3=infile.sym_p3
  p4=infile.sym_p4
  
  t=t1-t0
  if t gt 5400. then t=t-5400.*fix(t/5400)

if t lt 2700 then sm=0 else sm=p1*t+p2+t^2+p3*t^3+p4 
im_out=GAUSS_SMOOTH(im_in, sm, /EDGE_TRUNCATE)
return,im_out
end