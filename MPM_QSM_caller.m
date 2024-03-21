%%% Description: MPM QSM pipeline
% main steps:
% 1) complex-fit over echoes for pdw and t1w images,
%    simple phase difference for mtw images
%    for odd and even echoes done separately
% 2) ROMEO phase unwrapping
% 3) masking based on ROMEO quality map
% 4) rotation to scanner space
% 5) PDF background field removal
% 6) star QSM for dipole inversion as default (optional: non-linear dipole inversion)


%%% Publications:
% Please remember to give credit to the authors of the methods used:
% 1. SEPIA toolbox:
% Chan, K.-S., Marques, J.P., 2021. Neuroimage 227, 117611.
% 2. SPM12 - rigid body registration:
% Friston KJ, et al. Magnetic Resonance in Medicine 35 (1995):346-355
% 3. complex fit of the phase:
% Liu, Tian, et al. MRM 69.2 (2013): 467-476.
% 4. ROMEO phase uwnrapping:
% Dymerska, Barbara, and Eckstein, Korbinian et al. Magnetic Resonance in Medicine (2020).
% 5. PDF background field removal:
% Liu, Tian, et al. NMR in Biomed. 24.9 (2011): 1129-1136.
% 6. starQSM:
% Wei, Hongjiang, et al. NMR in Biomed. 28.10 (2015): 1294-1303.

%%% Inputs:
% romeo_command          : path to romeo phase uwnrapping followed by romeo command, i.e. (in linux) '/your_path/bin/romeo' or (in windows) 'D:\your_path\bin\romeo'
% in_root_dir            : root directory to input nifti files
% out_root_dir           : root directory to output nifti files
% B0                     : magnetic field strength, in Tesla
% dipole_inv             : dipole inversion method, either 'Star-QSM' or 'ndi'
%                          'ndi'      - non-linear dipole inversion
%                                       (also known as iterative Tikhonov),
%                                       may give more contrast than Star-QSM but is less robust to noise
%                          'Star-QSM' - is very robust to noise and quick
% calc_mean_qsm          : 'yes' or 'no', if 'yes' mean QSM from the three
%                           MPM acquisitions will be calculated, ATTENTION: currently nor
%                           coregistration between the scans is implemented

%%%% Inputs - directories, parameters and files specific to given contrast
% ATTENTION: ensure only niftis you want to use are in that folder, with increasing echo numbering:
% mag_dir                : % folder with magnitude niftis
% ph_dir                 : % folder with phase inftis
% TEs                    : % echo time in ms
% output_dir             : % output QSM directory for a specific MPM contrast
% mask_method            : % 'spm', 'spm_quality', 'quality', or 'load'        
% spm_mask_thr           : mask threshold between 0 and 1 used in 'spm', 'spm_quality'
% qmask_thr              : mask threshold between 0 and 1 used in 'spm_quality', 'quality'
% mask_file              : mask file used if mask_method = 'load'
% calc_mean_qsm          : % 'yes' or 'no' , if 'yes' it calculates mean QSM from all contrasts

%%% Outputs:
%%%% combined final results in out_root_dir:
% QSM_all_mean.nii             : mean QSM over all contrasts in scanner space (3rd dimension is along B0-axis)
% QSM_all_invrot_mean.nii      : mean QSM over all contrasts in image space (as acquired, for comparison with MPM quantitative maps)
% QSM_pdw_t1w_mean.nii         : mean QSM over PDw and T1w contrasts (without noisy MTw) in scanner space
% QSM_pdw_t1w_invrot_mean.nii  : mean QSM over PDw and T1w contrasts in image space

%%%% final results - per contrast in subfolders in out_root_dir:
% sepia_QSM.nii             : QSM in scanner space
% sepia_QSM_invrot.nii      : QSM in image space

%%%% additional outputs:
% ph.nii                    : two volumes (odd and even) of fitted phase
% ph_romeo.nii              : ph.nii unwrapped with ROMEO
% quality.nii               : quality map calculated by ROMEO algorithm and used for masking
% mask.nii                  : binary mask in image space
% mask_rot.nii              : binary mask in scanner space
% B0.nii                    : field map in Hz in image space
% B0_rot.nii                : field map in Hz in scanner space
% sepia_local-field.nii.gz  : map of local field variations (after background field removal using PDF)
% settings_romeo.txt        : settings used for ROMEO unwrapping (useful if unwrapping again outside MPM QSM the pipeline)
% header_sepia.mat          : header used for SEPIA toolbox (useful when exploring SEPIA GUI)

% script created by Barbara Dymerska
% @ UCL FIL Physics
% last modifications 09/09/2021

totstart = tic ;

%%%%% USER PARAMETERS %%%%%
para.romeo_command = 'romeo' ;
para.in_root_dir = '/home/bdymerska/Documents/data/3T/Gabor_Perlaki/pdw_QSM_data' ;
para.out_root_dir = '/home/bdymerska/Documents/data/3T/Gabor_Perlaki/pdw_QSM_data/';

para.mask_method= 'spm_quality';
para.spm_mask_thr=0.5;
para.qmask_thr=0.5;
para.mask_file = 'you/can/load/your/mask/if/you/want.nii';
para.B0 = 3;
para.dipole_inv = 'Star-QSM' ;




% directories, parameters and files specific to given contrast 
% ensure they are in the right order (PDw, T1w, MTw) 
% otherwise mean PDw+T1w QSM will be something different when calc_mean_qsm = 'yes':
for run = 1
    
    switch run
        case 1 %PDw
            para.mag_dir = 'pdw_kp_mtflash3d_v1s_caipi_0010' ; % folder with magnitude niftis
            para.ph_dir = 'pdw_kp_mtflash3d_v1s_caipi_RR_0022' ; % folder with phase inftis
            para.TEs = 2.3*[1 2 3 4 5 6 7 8] ; % echo time in ms
            para.output_dir = 'QSM_pdw' ; % output QSM directory for a specific MPM contrast
            
        case 2 % T1w
            para.mag_dir = 't1w_mfc_3dflash_v1k_0022' ;
            para.ph_dir = 't1w_mfc_3dflash_v1k_0023' ;
            para.TEs = [2.3 4.68 7.06 9.44 11.82 14.2] ;
            para.output_dir = 't1w_22_23' ;
            
        case 3 % MTw
            para.mag_dir = 'mtw_mfc_3dflash_v1k_180deg_0031' ;
            para.ph_dir = 'mtw_mfc_3dflash_v1k_180deg_0032' ;
            para.TEs = [2.2 4.58 6.96 9.34] ; 
            para.output_dir = 'mtw_31_32' ;
 
    end
    %%%%% END OF USER PARAMETERS %%%%%
    
    
    [QSM_V, QSM_all(:,:,:,run), QSMinvrot_V, QSM_all_invrot(:,:,:,run)] = MPM_QSM(para) ;
    
end


sprintf('total processing finished after %s' , secs2hms(toc(totstart)))