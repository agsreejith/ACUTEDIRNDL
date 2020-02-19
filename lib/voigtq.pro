function voigtq,wave, absr, line

  bnorm = absr.b/299792.4581
  vd = absr.b/ (line.wave * 1.0e-13)

  vel = abs((wave/(line.wave*(1.0+absr.z)) - 1.0)  / bnorm)
  a = line.gamma / (12.56637 * vd)

  calc1 = where(vel GE 10.0)
  calc2 = where(vel LT 10.0)

  vo = vel*0.0
  IF (calc1[0] NE -1) then begin
    vel2 = vel[calc1]*vel[calc1]
    hh1 = 0.56419/vel2 + 0.846/(vel2*vel2)
    hh3 = -0.56 / (vel2 * vel2)
    vo[calc1] = a * (hh1 + a * a * hh3)
  ENDIF
  if (calc2[0] NE -1) then vo[calc2] = voigt(a,vel[calc2])

  tau = 0.014971475*(10.0^absr.n)*line.f*vo/vd

  return, exp(-tau)
end
