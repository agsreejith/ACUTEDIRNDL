pro linearity,value,y

limit=60000
;f(x) = p1*x^3 + p2*x^2 + p3*x + p4
;Coefficients (with 95% confidence bounds):
;p1 =     0.03137  (-0.02162, 0.08437)
;p2 =     -0.4327  (-0.642, -0.2233)
;p3 =       1.475  (1.239, 1.712)
;p4 =     -0.1937  (-0.2662, -0.1212)

p1 = 0.03137
p2 = -0.4327
p3 = 1.475
p4 = -0.1937

if value gt 60000 then begin
  if value gt 72000 then y = 72000
  x=value/60000
  y=p1*x^3 + p2*x^2 + p3*x + p4
  y=y+60000
endif else begin
  y=value
endelse
end