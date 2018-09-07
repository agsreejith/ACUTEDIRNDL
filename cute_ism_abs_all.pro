function cute_ism_abs_all,flux,n_mg2,n_mg1,n_fe2
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


;;Parameters for MgI
;TCa2=7.7d5        ;temperature of the gas at the formation of the Ca2 emission in K
;mCa2=40.078/N_A  ;atomic mass of Ca2 in g: atomic mass in g/mole / Avogrado number

MgIw      = 2852.127 
MgI_loggf = 0.255
MgI_stark = -5.640

;;Parameters for FeII
;TCa2=7.7d5        ;temperature of the gas at the formation of the Ca2 emission in K
;mCa2=40.078/N_A  ;atomic mass of Ca2 in g: atomic mass in g/mole / Avogrado number

FeIIw      = 2599.39515
FeII_loggf = 0.378
FeII_stark = -6.53

;===============================
;ISM fixed parameters
ISM_b_Mg2=2.0        ;b-parameter for the Ca2 ISM lines in km/s
vr_ISM=0.          ;radial velocity of the ISM absorption lines in km/s

;===============================
;Construct the ISM absorption

n_flux=flux[1,*]/flux[2,*]

;for MgII doublet
absorberMg1=create_struct('ion','MG21','N',n_mg2,'B',ISM_b_Mg2,'Z',0.0)
lineMg1=create_struct('ion','Mg21','wave',MgII1w+MgII1w*vr_ISM/vc,'F',10^MgII1_loggf,'gamma',10^MgII1_stark)
ISMMg21=voigtq(flux[0,*],absorberMg1,lineMg1)

absorberMg2=create_struct('ion','MG22','N',n_mg2,'B',ISM_b_Mg2,'Z',0.0)
lineMg2=create_struct('ion','Mg22','wave',MgII2w+MgII2w*vr_ISM/vc,'F',10^MgII2_loggf,'gamma',10^MgII2_stark)
ISMMg22=voigtq(flux[0,*],absorberMg2,lineMg2)

;for MgI
absorberMgI=create_struct('ion','MG1','N',n_mg1,'B',ISM_b_Mg2,'Z',0.0)
lineMgI=create_struct('ion','Mg1','wave',MgIw*vr_ISM/vc,'F',10^MgI_loggf,'gamma',10^MgI_stark)
ISMMg1=voigtq(flux[0,*],absorberMgI,lineMgI)

;for FeII 
absorberFeII=create_struct('ion','FE"','N',n_fe2,'B',ISM_b_Mg2,'Z',0.0)
lineFeII=create_struct('ion','Mg1','wave',FeIIw*vr_ISM/vc,'F',10^FeII_loggf,'gamma',10^FeII_stark)
ISMFe2=voigtq(flux[0,*],absorberFeII,lineFeII)

;absorberCaK=create_struct('ion','Ca2K','N',N,'B',ISM_b_Ca2,'Z',0.0)
;lineCaK=create_struct('ion','Ca2K','wave',CaKwl0+CaKwl0*vr_ISM/vc,'F',10^CaK_loggf,'gamma',10^CaK_stark)
;ISMCaK=voigtq(flux[0,*],absorberCaK,lineCaK)
;
;absorberCaH=create_struct('ion','Ca2H','N',N,'B',ISM_b_Ca2,'Z',0.0)
;lineCaH=create_struct('ion','Ca2H','wave',CaHwl0+CaHwl0*vr_ISM/vc,'F',10^CaH_loggf,'gamma',10^CaH_stark)
;ISMCaH=voigtq(flux[0,*],absorberCaH,lineCaH)

ISM = ISMMg21*ISMMg22*ISMMg1*ISMFe2

flux_absorption = ISM * n_flux
flux[1,*]=flux[2,*]*flux_absorption

return,flux

end