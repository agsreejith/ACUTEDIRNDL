; NAME:
;      CUTE_LCA
;
; PURPOSE:
;      Inputs line core emission and ISM absorption in the stellar spectrum for MgII.
; CALLING SEQUENCE:
;      Result = COSMIC(flux,sigmaMg22,sigmaMg21,N,E)
;
; INPUTS:
;      flux       = The input array containing 3 columns, first column is wavleength, second is flux and third is continium
;      SigmaMg22  = The ISM width for second MgII line. 
;      SigmaMg21  = The ISM width for second MgII line.
;      N          = The Column density of MgII line.
;      E          = The stellar chromospheric emission at 1 AU.
; OUTPUT:
;      Returns the flux values with line core emission and ISM absorption added.
; REQUIRES:
;      The input flux, width of lines, column density and stellar chromospheric emission at 1 AU. 
;
; PROCEDURE:
;      Injects ine core emission and ISM absorption in the stellar spectrum for MgII. 
; MODIFICATION HISTORY:
;      created 06.12.2017 by A. G. Sreejith
;
;==============================================================================================
function voigtq,wave, abs, line

  bnorm = abs.b/299792.4581
  vd = abs.b/ (line.wave * 1.0e-13)

  vel = abs((wave/(line.wave*(1.0+abs.z)) - 1.0)  / bnorm)
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

  tau = 0.014971475*(10.0^abs.n)*line.f*vo/vd
  return, exp(-tau)
end
;==============================================================================================
function cute_lca,flux,sigmaMg22,sigmaMg21,N,E

;Constants
R_sun=6.957d10          ; cm
AU=1.49598d13           ; cm (Astronomical Unit)
vc=299792.458d0       ;speed of light in km/s
c0=2.99792458d10        ; cm/s (speed of light)
sigma=5.67051d-5        ; erg/cm^2/s/K^4 (stefan-boltzmann)
k_B=1.380658d-16        ; erg/K = g cm^2/s^2/K (boltzmann const)
N_A=6.02214179d23       ; /mol (Avagadro constant)
cA=2.99792458d18        ; Angstrom/s (speed of light)

;===============================

;===============================

;;Parameters for MgII
;TCa2=7.7d5        ;temperature of the gas at the formation of the Ca2 emission in K
;mCa2=40.078/N_A  ;atomic mass of Ca2 in g: atomic mass in g/mole / Avogrado number

MgII1w      = 2795.5280
MgII1_loggf = 0.100
MgII1_stark = -5.680    

MgII2w      = 2802.7050
MgII2_loggf = -0.210
MgII2_stark =-5.680     ;Stark damping constant

Mgaratio_loggf2to1=(10^MgII2_loggf)/(10^MgII1_loggf)
;;===============================
in_flux=flux[1,*]/flux[2,*]
fl=flux[1,*]
plot,flux[0,*],flux[1,*],xrange=[2595,3203]
WRITE_PNG, '1.png', TVRD(/TRUE)
;;===============================
;CaHwl0=3968.470    ;CaH wavelength
;CaH_loggf=-0.2      ;CaH loggf
;CaH_stark=8.193     ;Stark damping constant
;
;CaKwl0=3933.664    ;CaK wavelength
;CaK_loggf=0.105     ;CaK loggf
;CaK_stark=8.207     ;Stark damping constant
;
;Caratio_loggfKtoH=10^CaK_loggf/10^CaH_loggf
;;===============================
;Construct the Gaussian emission and measure the S and logR' values without ISM absorption
;It is not possible to use the S and logR' for this because the S caluclated in this way is not in the MW system
E=E*AU^2

Mg21em=E/(1.+Mgaratio_loggf2to1)
Mg22em=Mgaratio_loggf2to1*E/(1.+Mgaratio_loggf2to1)

gaussMg22=gaussian(flux[0,*],[0.3989*Mg22em/sigmaMg22,MgII2w,sigmaMg22])

gaussMg21=gaussian(flux[0,*],[0.3989*Mg21em/sigmaMg21,MgII1w,sigmaMg21])



;Ca2em=E/(1.+Caratio_loggfKtoH)
;Ca2Kem=Caratio_loggfKtoH*E/(1.+Caratio_loggfKtoH)
;
;sigmaCa2K=sqrt(k_B*TCa2/mCa2/c0^2)*CaKwl0; & sigmaCa2K=8.*sigmaCa2K
;gaussCa2K=gaussian(flux[0,*],[0.3989*Ca2Kem/sigmaCa2K,CaKwl0,sigmaCa2K])
;
;sigmaCa2H=sqrt(k_B*TCa2/mCa2/c0^2)*CaHwl0; & sigmaCa2H=8.*sigmaCa2H
;gaussCa2H=gaussian(flux[0,*],[0.3989*Ca2Hem/sigmaCa2H,CaHwl0,sigmaCa2H])

gaussMg2 = gaussMg21 + gaussMg22

flux_emission = flux[1,*] + gaussMg2
;===============================
;ISM fixed parameters
ISM_b_Mg2=2.0        ;b-parameter for the Ca2 ISM lines in km/s
vr_ISM=0.          ;radial velocity of the ISM absorption lines in km/s

;===============================
;Construct the ISM absorption 

n_flux=flux_emission/flux[2,*]

oplot,flux[0,*],flux_emission,psym=1
WRITE_PNG, '2.png', TVRD(/TRUE)

absorberMg1=create_struct('ion','MG21','N',N,'B',ISM_b_Mg2,'Z',0.0)
lineMg1=create_struct('ion','Mg21','wave',MgII1w+MgII1w*vr_ISM/vc,'F',10^MgII1_loggf,'gamma',10^MgII1_stark)
ISMMg21=voigtq(flux[0,*],absorberMg1,lineMg1)

absorberMg2=create_struct('ion','MG22','N',N,'B',ISM_b_Mg2,'Z',0.0)
lineMg2=create_struct('ion','Mg22','wave',MgII2w+MgII2w*vr_ISM/vc,'F',10^MgII2_loggf,'gamma',10^MgII2_stark)
ISMMg22=voigtq(flux[0,*],absorberMg2,lineMg2)

;absorberCaK=create_struct('ion','Ca2K','N',N,'B',ISM_b_Ca2,'Z',0.0)
;lineCaK=create_struct('ion','Ca2K','wave',CaKwl0+CaKwl0*vr_ISM/vc,'F',10^CaK_loggf,'gamma',10^CaK_stark)
;ISMCaK=voigtq(flux[0,*],absorberCaK,lineCaK)
;
;absorberCaH=create_struct('ion','Ca2H','N',N,'B',ISM_b_Ca2,'Z',0.0)
;lineCaH=create_struct('ion','Ca2H','wave',CaHwl0+CaHwl0*vr_ISM/vc,'F',10^CaH_loggf,'gamma',10^CaH_stark)
;ISMCaH=voigtq(flux[0,*],absorberCaH,lineCaH)

ISM = ISMMg21*ISMMg22

flux_absorption = ISM * n_flux


flux[1,*]=flux[2,*]*flux_absorption

writecol,'testmg.txt',flux[0,*],fl,flux_emission,flux[1,*]
oplot,flux[0,*],flux[1,*],psym=5
WRITE_PNG, '3.png', TVRD(/TRUE)
;===============================
return,flux
end
