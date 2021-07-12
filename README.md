# MPM_QSM
QSM pipeline for Multi Parametric Mapping acquisitions


# Main computational steps:

 1) complex-fit over echoes for pdw and t1w images,
    simple phase difference for mtw images
    for odd and even echoes done separately
 2) ROMEO phase unwrapping
 3) masking based on ROMEO quality map
 4) rotation to scanner space for oblique acquisitions
 5) PDF background field removal
 6) star QSM for dipole inversion as default (optional: non-linear dipole inversion)


# Installation steps:

1. Download Zip with all the files from:

	https://github.com/baskadym/MPM_QSM.git
	or clone it to your GitHub repository

2. Download compiled version of ROMEO either for windows or linux:
	
	https://github.com/korbinian90/ROMEO/releases/tag/v3.2.0
	and unzip it in chosen destination

3. Download sepia toolbox:
	
	https://github.com/kschan0214/sepia.git
	and unzip it in chosen destination

4. Set you local paths to MEDI and STI toolboxes downloaded in step 1 in file:
 	/your_path/sepia-master/SpecifyToolboxesDirectory.m

	as following:
	MEDI_dir = '/your/MEDI/path/';
	STISuite_dir = '/your/STI/path/';
	FANSI_dir = [];
	SEGUE_dir = [];

5. Add to your matlab path: SEPIA toolbox and MPM_QSM folder
6. Edit MPM_QSM.m user parameters, where you specify folders to you nifti files


# Publications:

Please remember to give credit to the authors of the methods used:

1. SEPIA toolbox:
Chan, K.-S., Marques, J.P., 2021. SEPIA—Susceptibility mapping pipeline tool for phase images. Neuroimage 227, 117611.

2. complex fit of the phase:
Liu, Tian, et al. "Nonlinear formulation of the magnetic field to source relationship for robust quantitative susceptibility mapping." Magnetic resonance in medicine 69.2 (2013): 467-476.

3. ROMEO phase uwnrapping:
Dymerska, Barbara, and Eckstein, Korbinian et al. "Phase unwrapping with a rapid opensource minimum spanning tree algorithm (ROMEO)." Magnetic Resonance in Medicine (2020).

4. PDF background field removal:
Liu, Tian, et al. "A novel background field removal method for MRI using projection onto dipole fields." NMR in Biomedicine 24.9 (2011): 1129-1136.

5. starQSM:
Wei, Hongjiang, et al. "Streaking artifact reduction for quantitative susceptibility mapping of sources with large dynamic range." NMR in Biomedicine 28.10 (2015): 1294-1303.

