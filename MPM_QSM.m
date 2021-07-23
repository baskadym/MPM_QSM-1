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
% 2. complex fit of the phase:
% Liu, Tian, et al. MRM 69.2 (2013): 467-476.
% 3. ROMEO phase uwnrapping:
% Dymerska, Barbara, and Eckstein, Korbinian et al. Magnetic Resonance in Medicine (2020).
% 4. PDF background field removal:
% Liu, Tian, et al. NMR in Biomed. 24.9 (2011): 1129-1136.
% 5. starQSM:
% Wei, Hongjiang, et al. NMR in Biomed. 28.10 (2015): 1294-1303.

%%% Inputs:
% romeo_command          : path to romeo phase uwnrapping followed by romeo command, i.e. (in linux) '/your_path/bin/romeo' or (in windows) 'D:\your_path\bin\romeo'
% B0                     : magnetic field strength, in Tesla
% algorParam.qsm.method  : dipole inversion method, either 'Star-QSM' or 'ndi'
%                          'ndi' - non-linear dipole inversion (also known as iterative Tikhonov), may give more contrast than Star-QSM but is less robust to noise
%                          'Star-QSM' is very robust to noise and quick
% in_root_dir            : root directory to input nifti files
% out_root_dir           : root directory to output nifti files
%%%% Inputs - directories, parameters and files specific to given contrast
% mag_dir                : % folder with magnitude niftis
% ph_dir                 : % folder with phase inftis
% TEs                    : % echo time in ms
% output_dir             : % output directory for a specific submeasurement from MPM
% mag_file               : % magnitude reference nifti file for ROMEO unwrapping and masking


%%% Outputs:
%%%% combined final results in out_root_dir:
% QSM_all_mean.nii             : mean QSM over all contrasts in scanner space (3rd dimension is along B0-axis)
% QSM_all_invrot_mean.nii      : mean QSM over all contrasts in image space (as acquired, for comparison with MPM quantitative maps)
% QSM_pdw_t1w_mean.nii         : mean QSM over PDw and T1w contrasts (without noisy MTw) in scanner space
% QSM_pdw_t1w_invrot_mean.nii  : mean QSM over PDw and T1w contrasts in image space

%%%% final results - per contrast in subfolders in out_root_dir:
% sepia_QSM.nii.gz          : QSM in scanner space
% sepia_QSM_invrot.nii.gz   : QSM in image space

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
% last modifications 12/07/2021

tstart = tic ;
%%%%% USER PARAMETERS %%%%%
romeo_command = '~/Documents/MRI_software/ROMEO/romeo_linux_3.2.0/bin/romeo' ;

B0 = 7;
algorParam.qsm.method = 'Star-QSM' ;

in_root_dir = '/media/barbara/hdd2/DATA/FIL/7T/20210623.M700198_FIL_analysis' ;
out_root_dir = '/media/barbara/hdd2/DATA/FIL/7T/20210623.M700198_FIL_analysis/SEPIA/MORSE_scan2';


% directories, parameters and files specific to given contrast:
for run = 1:3
    
    switch run
        case 1 %pdw
            mag_dir = 'pdw_mfc_3dflash_v1k_0025' ; % folder with magnitude niftis
            ph_dir = 'pdw_mfc_3dflash_v1k_0026' ; % folder with phase inftis
            TEs = [2.2 4.58 6.96 9.34 11.72 14.1] ; % echo time in ms
            output_dir = 'pdw_25_26' ; % output directory for a specific submeasurement from MPM
            mag_file = 's2021-06-23_10-18-112654-00001-01728-6.nii' ; % magnitude reference nifti file for ROMEO unwrapping and masking
            
        case 2 % t1w
            mag_dir = 't1w_mfc_3dflash_v1k_0022' ;
            ph_dir = 't1w_mfc_3dflash_v1k_0023' ;
            TEs = [2.3 4.68 7.06 9.44 11.82 14.2] ;
            output_dir = 't1w_22_23' ;
            mag_file = 's2021-06-23_10-18-111631-00001-01728-6.nii' ;
            
        case 3 % mtw
            mag_dir = 'mtw_mfc_3dflash_v1k_180deg_0031' ;
            ph_dir = 'mtw_mfc_3dflash_v1k_180deg_0032' ;
            TEs = [2.2 4.58 6.96 9.34] ; % echo time in ms
            output_dir = 'mtw_31_32' ;
            mag_file = 's2021-06-23_10-18-114211-00001-01152-4.nii' ; % magnitude reference nifti file for ROMEO unwrapping and masking
    end
    
    
    %%%%% END OF USER PARAMETERS %%%%%
    
    mag_fulldir = fullfile(in_root_dir, mag_dir) ;
    ph_fulldir = fullfile(in_root_dir, ph_dir) ;
    
    output_fulldir = fullfile(out_root_dir, output_dir) ;
    if ~exist(output_fulldir, 'dir')
        mkdir(output_fulldir)
    end
    cd(output_fulldir)
    
    TEs = TEs/10^3 ;
    
    ph_files = dir(ph_fulldir);
    mag_files = dir(mag_fulldir);
    
    for t = 1:size(TEs,2)
        
        ph_1tp = load_untouch_nii(fullfile(ph_fulldir, ph_files(t+2).name));
        ph(:,:,:,t) = ph_1tp.img ;
        
        if size(TEs,2) >= 6 % for mtw acquisition we cannot perform complex fit so we don't need magnitude data for each TE
            
            mag_1tp = load_untouch_nii(fullfile(mag_fulldir, mag_files(t+2).name)) ;
            mag(:,:,:,t) = mag_1tp.img ;
            
        end
        
    end
    ph_1tp_spm = nifti(fullfile(ph_fulldir, ph_files(t+2).name)) ;
    Maff_orig = ph_1tp_spm.mat ;
    O = Maff_orig\[0 0 0 1]' ;
    O = O(1:3)' ;
    
    clear ph_1tp.img mag_1tp
    
    % rescaling phase into [0,2*pi] phase range
    ph = 2*pi*single(ph - min(vector(ph)))/single(max(vector(ph))-min(vector(ph))) ;
    
    for read_dir = 1:2
        
        if size(TEs,2) < 6 % complex fit is only possible if at least 3 echoes per readout direction available
            disp('calculating phase difference')
            FM = angle(exp(1i*(ph(:,:,:,read_dir+2)-ph(:,:,:,read_dir)))) ;
        else
            disp('complex fitting phase')
            compl = single(mag).*exp(-1i*ph);
            [FM, ~, ~, ~] = Fit_ppm_complex_TE(compl(:,:,:,read_dir:2:end),TEs(read_dir:2:end));
        end
        
        FM_both(:,:,:,read_dir) = FM ;
        
    end
    clear mag ph FM
    
    % saving odd and even echoes as one file for ROMEO unwrapping
    FM_file = 'ph.nii' ;
    FM_both = make_nii(FM_both) ;
    FM_both.hdr.hist = ph_1tp.hdr.hist ;
    
    centre_and_save_nii(FM_both, FM_file, ph_1tp.hdr.dime.pixdim);
    clear FM_both
    
    % creating magnitude reference for ROMEO unwrapping
    mag_fullfile = fullfile(mag_fulldir, mag_file);
    mag_ref = load_untouch_nii(mag_fullfile) ;
    mag_ref = repmat(mag_ref.img,[1 1 1 2]) ;
    mag_fullfile = sprintf('mag_TE%i.nii',size(TEs,2)) ;
    centre_and_save_nii(make_nii(mag_ref), mag_fullfile, ph_1tp.hdr.dime.pixdim);
    clear mag_ref
    
    disp('phase unwrapping with ROMEO and:')
    disp('...removing global mean value')
    disp('......field map calculation')
    disp('.........saving quality map for masking')
    TE = (TEs(3)-TEs(1))*10^3; % effective echo time difference after phase complex fitting in seconds
    [~, FM_name,~] = fileparts(FM_file) ;
    FM_romeo_file = sprintf('%s_romeo.nii',FM_name) ;
    if isunix
        status =  unix(sprintf('%s %s -m %s -o %s -t [%i,%i] -k nomask -g -q -B', romeo_command, FM_file, mag_fullfile, FM_romeo_file, TE, TE)) ;
    elseif ispc
        status = system(sprintf('%s %s -m %s -o %s -t [%i,%i] -k nomask -g -q -B', romeo_command, FM_file, mag_fullfile, FM_romeo_file, TE, TE)) ;
    end
    if status == 1
        error('ROMEO did not run properly - check your installation path')
    end
    delete(sprintf('mag_TE%i.nii',size(TEs,2)));
    %%%%%%
    % defining affine matrix in scanner space for data rotation to scanner
    % space with mantaining the same image origin (i.e. no translation)
    pixdim = ph_1tp.hdr.dime.pixdim(2:4) ;
    data_dim = ph_1tp.hdr.dime.dim(2:4) ;
    M_scanner(1,:) = [0 0 pixdim(3) -pixdim(3)*O(3)] ;
    M_scanner(2,:) = [pixdim(1) 0 0 -pixdim(1)*O(1)] ;
    M_scanner(3,:) = [0 pixdim(2) 0 -pixdim(2)*O(2)] ;
    M_scanner(4,:) = [0 0 0 1] ;
    
    FM = nifti('B0.nii') ;
    fm_data = FM.dat(:,:,:) ;
    fm_data(isnan(fm_data)) = 0 ;
    FM_V = spm_vol('B0.nii');
    spm_write_vol(FM_V, fm_data) ;
    
    img2scanner_mat = FM.mat\M_scanner ;
    FMrot = zeros(data_dim) ;
    
    for slice = 1 : data_dim(3)
        FMrot(:,:,slice) = spm_slice_vol(FM_V, img2scanner_mat*spm_matrix([0 0 slice]), data_dim(1:2), -7) ;
    end
    
    FMrot(isnan(FMrot)) = 0 ;
    FMrot_V = FM_V ;
    FMrot_V.mat = M_scanner ;
    FMrot_V.fname = 'B0_rot.nii';
    spm_write_vol(FMrot_V, FMrot)

    
    disp('quality masking')
    qmap = load_untouch_nii('quality.nii') ;
    qmap_bin = qmap.img ;
    qmap_bin(qmap.img>0.3) = 1 ;
    qmap_bin(qmap.img<=0.3) = 0 ;
    clear qmap
    qmap_bin(isnan(qmap_bin)) = 0 ;
    qmap_bin = imfill(qmap_bin,6,'holes') ;
    qmask = smooth3(qmap_bin, 'gaussian') ;
    qmask(qmask>0.6) = 1 ;
    qmask(qmask<=0.6) = 0 ;
    clear qmap_bin
    
    qmask = make_nii(int16(qmask)) ;
    qmask.hdr.hist = ph_1tp.hdr.hist ;
    qmask.hdr.dime.pixdim = ph_1tp.hdr.dime.pixdim ;
    qmask_file = fullfile(output_fulldir, 'mask.nii') ;
    save_nii(qmask, qmask_file);
    
    mask_rot = zeros(data_dim) ;
    mask_V = spm_vol('mask.nii');
    for slice = 1 : data_dim(3)
        mask_rot(:,:,slice) = spm_slice_vol(mask_V, img2scanner_mat*spm_matrix([0 0 slice]), data_dim(1:2), -7) ;
    end
    mask_rot(isnan(mask_rot)) = 0 ;
    mask_rot_V = mask_V ;
    mask_rot_V.mat = M_scanner ;
    mask_rot_V.fname = 'mask_rot.nii';
    spm_write_vol(mask_rot_V, mask_rot)

    clear qmask
    
    %% SEPIA - calculates QSM
    
    % create SEPIA header
    CF = B0*42.58*1e6;	% imaging frequency, in Hz (B0*gyromagnetic_ratio*1e6)
    delta_TE = 1;	    % echo spacing, in second - we have already combined data, in such situation set to 1
    B0_dir = [0;1;0];	% main magnetic field direction, it's always [0,1,0] because the images are resliced so that 2nd dimension is aligned with B0
    hdr = load_nii_hdr('B0_rot.nii') ;
    matrixSize = hdr.dime.dim(2:4) ;	    % image matrix size
    voxelSize = ph_1tp.hdr.dime.pixdim(2:4) ;	% spatial resolution of the data, in mm
    header_fullfile = fullfile(output_fulldir, 'header_sepia.mat') ;
    save(header_fullfile, 'B0', 'B0_dir', 'CF', 'TE', 'delta_TE', 'matrixSize', 'voxelSize')
    
    % general SEPIA parameters
    sepia_addpath
    
    algorParam.general.isBET = 0 ;
    algorParam.general.isInvert = 1 ;
    algorParam.general.isGPU = 0 ;
    
    % inputs for background field removal
    input(1).name = 'B0_rot.nii' ;
    input(2).name = 'mask_rot.nii' ;
    input(4).name = header_fullfile ;
    
    algorParam.bfr.refine = 0 ;
    algorParam.bfr.erode_radius = 0 ;
    algorParam.bfr.method = 'pdf' ;
    algorParam.bfr.tol = 0.1 ;
    algorParam.bfr.iteration = 50 ;
    algorParam.bfr.padSize = 30 ;
    
    
    % inputs for dipole inversion
    if strcmp(algorParam.qsm.method , 'ndi')
        
        algorParam.qsm.method = 'ndi' ;
        algorParam.qsm.tol = 1 ;
        algorParam.qsm.maxiter = 200 ;
        algorParam.qsm.stepSize = 1 ;
        
    elseif strcmp(algorParam.qsm.method , 'Star-QSM')
        
        algorParam.qsm.padsize = ones(1,3)*12 ;
        
    end
    
    output_basename = fullfile(output_fulldir, 'sepia') ;
    
    disp('background field removal')
    BackgroundRemovalMacroIOWrapper(input,output_basename,input(2).name,algorParam);
    
    disp('dipole inversion')
    input(1).name = fullfile(output_fulldir, 'sepia_local-field.nii.gz') ;
    QSMMacroIOWrapper(input,output_basename,input(2).name,algorParam);
    
    sprintf('run %i finished after %s' ,run, secs2hms(toc(tstart)))
    
    QSM_V = spm_vol(fullfile(output_fulldir, 'sepia_QSM.nii.gz'));
    
    disp('rotation of QSM back to the original image space')
    
    %%%%%%
    gunzip('sepia_QSM.nii.gz') 
    QSM = nifti('sepia_QSM.nii') ;
    
    scanner2img_mat = M_scanner\FM.mat ;
    QSMinvrot = zeros(data_dim) ;
    QSM_V = spm_vol('sepia_QSM.nii');
    for slice = 1 : data_dim(3)
        QSMinvrot(:,:,slice) = spm_slice_vol(QSM_V, scanner2img_mat*spm_matrix([0 0 slice]), data_dim(1:2), -7) ;
    end
    
    QSMinvrot_V = FM_V ;
    QSMinvrot_V.fname = 'sepia_QSM_invrot.nii';
    spm_write_vol(QSMinvrot_V, QSMinvrot)
 
    QSM_all(:,:,:,run) = QSM.dat(:,:,:) ;
    QSM_all_invrot(:,:,:,run) = QSMinvrot ;
    
    %     clear QSM QSM_invrot
    delete('sepia_mask-qsm.nii.gz')
    
end

QSM_all_mean = mean(QSM_all, 4) ;
QSM_pdw_t1w_mean = mean(QSM_all(:,:,:,1:2), 4) ;

QSM_all_invrot_mean = mean(QSM_all_invrot, 4) ;
QSM_pdw_t1w_invrot_mean = mean(QSM_all_invrot(:,:,:,1:2), 4) ;

QSM_all_invrot_mean = make_nii(QSM_all_invrot_mean, ph_1tp.hdr.dime.pixdim(2:4)) ;
% QSM_all_invrot_mean.hdr.hist = ph_1tp.hdr.hist ;

QSM_pdw_t1w_invrot_mean = make_nii(QSM_pdw_t1w_invrot_mean, ph_1tp.hdr.dime.pixdim(2:4)) ;
% QSM_pdw_t1w_invrot_mean.hdr.hist = ph_1tp.hdr.hist ;

centre_and_save_nii(make_nii(QSM_all_mean), fullfile(out_root_dir,'QSM_all_mean.nii'), ph_1tp.hdr.dime.pixdim);
centre_and_save_nii(make_nii(QSM_pdw_t1w_mean), fullfile(out_root_dir,'QSM_pdw_t1w_mean.nii'), ph_1tp.hdr.dime.pixdim);

centre_and_save_nii(QSM_all_invrot_mean, fullfile(out_root_dir,'QSM_all_invrot_mean.nii'), ph_1tp.hdr.dime.pixdim);
centre_and_save_nii(QSM_pdw_t1w_invrot_mean, fullfile(out_root_dir,'QSM_pdw_t1w_invrot_mean.nii'), ph_1tp.hdr.dime.pixdim);


sprintf('total processing finished after %s' , secs2hms(toc(tstart)))
