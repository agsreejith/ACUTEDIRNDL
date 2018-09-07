

README file for ACUTEDIRNDL
---------------------------

0. FEATURES
-----------

ACUTEDIRNDL is an IDL tool capable to create simulated data for the CUTE spectrograph, also considering the observation of NUV planetary transits.

Output will be in fits format
http://en.wikipedia.org/wiki/FITS


1. PRE-REQUISITIES
------------------

IDL astro packages.
gaussbroad.pro
gm_read_textstructure.pro
gm_read_textfile.pro
cmset_op.pro
get_tags.pro
keyword_defined.pro
mpfitfun

All required codes not present in the IDL astronomy library are provided in the library folder.

Add the simulation folder (and sub-folders) to the IDL source in the startup file.

2. USAGE
--------
This software is governed by an input parameter file called: cutedrndl_parameters.txt

The ACUTEDIRNDL routines are available in the src folder. Run cutedrndl.pro (just type cutedrndl inside 
the IDL environment) after editing the parameter file.
The outputs will be stored in the folder specified in input parameter file. The number of output files will 
depend on both transit duration and exposure time per observation. Note that the simulator overwrites the 
files in the output directory. The detailed description of the parameters in the input parameter file is described here below.

- File containing the stellar flux at earth. (This input can be commented out)
The stellar flux at Earth is to be provided in a file containing two columns, where the first is the wavelength
in Angstrom and the second in the flux in erg cm^-2 s^-1 A^-1. The file name has to contain the full path, 
thus directing the simulator to the location of the file.

- Stellar parameters (These inputs are required)
The stellar parameters required as input to the the simulator are 1) the stellar temperature (from 3500 to 
10000K, in steps of 100K); 2) the stellar radius (in unit of solar radii); 3) and Johnson V-band stellar_magnitude. 
It is also necessary to provide the location (i.e., path) of the files containing the stellar model fluxes 
that are provided together with the simulator. These files contain the model fluxes where the first column 
is the wavelength in Angstrom, the second column is the stellar flux in erg cm^-2 s^-1 Hz^-1, and the third 
column is the stellar continuum in erg cm^-2 s^-1 Hz^-1. Note that the stellar continuum is used for the 
inclusion of line core emission and interstellar medium absorption.

- Target observing parameters (These inputs are required)
The coordinates of the target in RA and DEC (in degrees) and exposure time per observation in seconds.	
The date of observation can be provided in the format of Julian date (JD). If JD is set to 0 then the time 
of the first observation will be taken from the system time, thus setting all following observations 
according to the exposure time and CCD readout time, that are both input parameters. If JD is not set to 
zero, then it should be set in Julian date and the code takes it as being the mid transit time. The time of 
all other frames is then computd accordingly following the transit parameters (see below), the exposure time,
 and CCD readout time.

- Instrument parameters (These inputs are required)
The instrument parameters include the position of the star on the slit and the spectrograph spectral 
resolution in Angstroms. The code allows for the position of the star on the slit from the slit center. The spectral 
resolution is a single value for the while wavelength range, but an update is planned in order to allow for 
a wavelength dependent spectral resolution. The slit dimaensions in arcseconds also have to be specified.

- Instrument specifications 1 (These inputs are required)
The first set of instrument specifications consists in the files containing the instrument wavelength scale,
 in Angstrom, and effective area, in cm^2. The first column is wavelength in anstroms and second column is effective area in cm^2.
 The user needs to specify also the full path of these files. For the CUTE spectrograph, these files are 
 provided in the extra directory. 

- Instrument specifications 2 (These inputs are required)
The second set of instrument specifications consists in the files containing the information on the slope of
 the spectrum on the CCD and the shape of the spectrum in the direction perpendicular to the dispertion 
 direction.The columns specify the shift in pixels from center and relative flux for different wavelength 
 groups. The user needs to specify also the full path of these files. For the CUTE spectrograph, these files
 are provided in the extra directory.

- Background parameters (These inputs are required)
The background parameters set how the background is calculated and whether to include background stars. The
 parameter "type_bckg" can be set equal to "default" or "calc". The first indicates that the background will
 be calculated in the same way as the "high bakcground" level set in the STIS data exposure time calculator. 
 he second indicates that the background will be computed on the basis of the actual zodiacal emission maps and the pointing 
 position of the satellite. To include background stars set "stars" equal to 1. In this case, 
 the simulator will query vizier to get the position of the stars around the region of the target and place them on their 
 respective position inside the slit.

- General stellar data (This input is required)
The stellar_params file, located in the "extra" directory, gives the main stellar parameters as a finction 
of spectral type. The first column is the spectral type, the second column is the stellar effective 
temperature, the third column is the B-V color, while the fourth column is the stellar radius in solar radii
. The user needs to specify also the full path of this file.

- Background file (This input is required)
The bckg_default file, located in the "extra" directory, gives gives the high background values in
STIS background interpolated and trimmed to CUTE wavelength band The user needs to specify also the full 
path of this file.

- Cosmic rays (This input is required)
This parameter sets the maximum number of cosmic rays present in each image. Each cosmic ray will create a 
sequence of 3 consecutive saturated pixels randomly placed and oriented across the CCD.

- Files for the zodiacal background calculations (This input is required)
Files indicating the parameters used for the zodiacal background calculations. The user needs to specify 
also the full path of these files.

- CCD parameters (This input is required)
These parameters specify the characteristics of the CCD. X_ and Y_pixels indicate the x and y pixel size of
 the CCD. The other parameters are the readout noise, in ADU, the dark level, in e/pixel/s , the average bia
 s level, in ADU, the readout time, in seconds [LUCA: SECONDS?], and the CCD gain, in photoelectrons per ADU . 

- Spacecraft parameters (This input is required)
The jitter_sim parameter sets whether spacecraft jitter has to be included or not (1 or 0, respectively). 
The jitter parameter sets the rms amount of spacecraft jitter in arcseconds. The rotation parameter indicates
the position angle (rotation with respect to the RA axis) of the slit, in degrees.
	
- Transit parameters (This input is required)
The impact parameter of the transit is between 0 and 1, where 0 corresponds to central transit.
 The planetary orbital period is in days, while the orbital inclination axis is in degrees. The parameter 
 rlt_semi_major_axis sets the orbital semi-major axis in units of stellar radii and mid transit time of 
 transit in JD.  The limb_file parameter sets the directory location of the PHOENIX models employed to 
 compute limb darkening coefficients as a function of stellar temperature and wavelength.

- Planetary radius  (One of these inputs is required)
The planetary radius has to be set either as a single value in units of stellar radii, in which case the
planetary radius will not depend on wavelength, or as a file containing wavelength (first column) and 
planetary radius as a function of wavelength (second column). Also in this case, the planetary radius as a 
function of wavelength has to be set in units of stellar radii. If the file with planetary radius as a 
functioon of wavelength is specified, then the code takes this instead of the constant radius.

Other parameters
- Currently the CCD linearity curve is set to the linear up to 60000 counts, after which it follows a 3rd 
order polynomial till 72000 counts (saturation). 
- The surface gravity and metallicity of the stars is assumed to be always equal to 4.5 and solar, 
respectively.	
	

3. DISCLAIMER
-------------
You may encounter bugs in this software. If you do, please report them. Your bug
reports are valuable contributions, since they allow us to notice and fix
problems on machines/platforms we don't have, and/or remained un-noticed.


4. REPORTING BUGS
-----------------
Drop in an email to: sreejith.aickara_AT_oeaw.ac.at

please include as many details as possible.

-----------------------------------------------------------
last modified: modified 07.09.2018 by A. G. Sreejith
