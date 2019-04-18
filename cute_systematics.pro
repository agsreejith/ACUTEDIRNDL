function cute_systematics,infile,t1,t0,im_in,type

;  f(x) = p1*x^3 + p2*x^2 + p3*x + p4
;  Coefficients (with 95% confidence bounds):
;  p1 =  -3.453e-12  (-6.928e-12, 2.124e-14)
;  p2 =   4.092e-08  (1.734e-10, 8.167e-08)
;  p3 =  -0.0001206  (-0.0002752, 3.391e-05)
;  p4 =         1.1  (0.9112, 1.289)
if keyword_defined(type) eq 0 then type = 'flux'

  s=size(im_in)
  nx=s[1]
  ny=s[2]
  im_out=dblarr(nx,ny)
  p1=infile.sym_p1
  p2=infile.sym_p2
  p3=infile.sym_p3
  p4=infile.sym_p4
  if infile.sc_orbit_period eq 0 then orbit = 5400. else orbit = double(infile.sc_orbit_period)*60
  t=(t1-t0)*86400 
  tprint=t
  if t ge orbit then t=t-orbit*fix(t/orbit)
  t=t/orbit
  sm=p1*t^3+p2*t^2+p3*t+p4
  logprint,tprint,sm,logfile=infile.file_out+'systematics.txt',logonly=1   
  if type eq 'defocus' then begin
    im_out=GAUSS_SMOOTH(im_in, sm, /EDGE_TRUNCATE)
    return,im_out
  endif else if type eq 'flux' then begin
    if sm eq 0.0 then im_out=im_in else im_out=sm*im_in
    return,im_out  
  endif
end