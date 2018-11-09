

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

	
3. EXAMPLES OF SOME INPUT FILES
-------------------------------	
Example for wavelength depended transit light curve

0	2000	2010	2020	2030	2040	2050
0	1	1	1	1	1	1

1	1	1	1	1	1	1

2	1	1	1	1	1	1

3	0.99	0.99	0.99	0.99	0.99	0.99

4	0.985	0.98	0.98	0.985	0.984	0.984

5	0.98	0.97	0.97	0.985	0.98	0.98

6	0.985	0.98	0.98	0.985	0.984	0.984

7	0.99	0.99	0.99	0.99	0.99	0.99

8	1	1	1	1	1	1

9	1	1	1	1	1	1

10	1	1	1	1	1	1



Example for input for shape of the footprint as a function of wavelength (first column is the pixels from centroid and other column is spread at different wavelength bins)

         -6.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0034 0.0345
         -5.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0001 0.0112 0.0330 0.0533
         -4.0000 0.0092 0.0048 0.0008 0.0002 0.0000 0.0015 0.0087 0.0335 0.0519 0.0606 0.1023
         -3.0000 0.0287 0.0372 0.0442 0.0298 0.0386 0.0499 0.0544 0.1544 0.2728 0.3066 0.3312
         -2.0000 0.0411 0.0541 0.0628 0.0514 0.0537 0.2941 0.4131 0.4800 0.4485 0.4074 0.4410
         -1.0000 0.0522 0.0680 0.1334 0.3568 0.5164 0.4159 0.3723 0.4148 0.4273 0.4088 0.4841
          0.0000 0.2114 0.4748 0.7110 0.3470 0.3495 0.3884 0.3997 0.5094 0.6004 0.6474 0.7763
          1.0000 0.4539 0.4948 0.4957 0.3768 0.4219 0.4596 0.4400 0.4548 0.3428 0.1823 0.0597
          0.0000 0.0244 0.0019 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000
          1.0000 0.0001 0.0004 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000       

4. DISCLAIMER
-------------
You may encounter bugs in this software. If you do, please report them. Your bug
reports are valuable contributions, since they allow us to notice and fix
problems on machines/platforms we don't have, and/or remained un-noticed.


5. REPORTING BUGS
-----------------
Drop in an email to: sreejith.aickara_AT_oeaw.ac.at

please include as many details as possible.

-----------------------------------------------------------
last modified: modified 09.11.2018 by A. G. Sreejith
