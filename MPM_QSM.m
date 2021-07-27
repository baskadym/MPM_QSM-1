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
% last modifications 27/07/2021

function [QSM , QSMinvrot] = MPM_QSM(para)

tstart = tic ;

    mag_fulldir = fullfile(para.in_root_dir, para.mag_dir) ;
    ph_fulldir = fullfile(para.in_root_dir, para.ph_dir) ;
    
    output_fulldir = fullfile(para.out_root_dir, para.output_dir) ;
    if ~exist(output_fulldir, 'dir')
        mkdir(output_fulldir)
    end
    cd(output_fulldir)
    
    TEs = para.TEs/10^3 ;
    
    ph_files = dir(ph_fulldir);
    mag_files = dir(mag_fulldir);
    
    for t = 1:size(TEs,2)
        
        ph_1tp = nifti(fullfile(ph_fulldir, ph_files(t+2).name));
        ph(:,:,:,t) = ph_1tp.dat(:,:,:) ;
        
        if size(TEs,2) >= 6 % for mtw acquisition we cannot perform complex fit so we don't need magnitude data for each TE
            
            mag_1tp = nifti(fullfile(mag_fulldir, mag_files(t+2).name)) ;
            mag(:,:,:,t) = mag_1tp.dat(:,:,:) ;
            
        end
        
    end
    
    
    % rescaling phase into [0,2*pi] phase range
    ph = 2*pi*single(ph - min(vector(ph)))/single(max(vector(ph))-min(vector(ph))) ;
    
    clear ph_V % otherwise if there is already this structure in memory it does not like it
    ph_V(1) = spm_vol(fullfile(ph_fulldir, ph_files(3).name)) ;
    ph_V(1).fname = 'ph.nii' ;
    ph_V(1).dt = [16 0] ;
    ph_V(1).descript = '';
    ph_V(2) =  ph_V(1) ;
    ph_V(2).n = [2 1] ;
    
    clear mag_V % otherwise if there is already this structure in memory it does not like it
    mag_V(1) = spm_vol(fullfile(mag_fulldir, mag_files(size(TEs,2)+2).name));
    mag_V(1).fname = sprintf('mag_TE%i.nii',size(TEs,2));
    mag_V(1).descript = '';
    mag_V(2) = mag_V(1) ;
    mag_V(2).n = [2 1] ;
    mag_ref = nifti(fullfile(mag_fulldir, mag_files(size(TEs,2)+2).name)) ;
    mag_ref = repmat(mag_ref.dat(:,:,:),[1 1 1 2]) ;
    for read_dir = 1:2
        
        if size(TEs,2) < 6 % complex fit is only possible if at least 3 echoes per readout direction available
            disp('calculating phase difference')
            FM(:,:,:,read_dir) = angle(exp(1i*(ph(:,:,:,read_dir+2)-ph(:,:,:,read_dir)))) ;
        else
            disp('complex fitting phase')
            compl = single(mag).*exp(-1i*ph);
            [FM_1, ~, ~, ~] = Fit_ppm_complex_TE(compl(:,:,:,read_dir:2:end),TEs(read_dir:2:end));
            FM(:,:,:,read_dir) = FM_1 ;
        end
        
        spm_write_vol(ph_V(read_dir),  FM(:,:,:,read_dir)) ;
        spm_write_vol(mag_V(read_dir), mag_ref(:,:,:,read_dir)) ;
        
    end
    clear mag ph FM mag_ref
    
    
    disp('phase unwrapping with ROMEO and:')
    disp('...removing global mean value')
    disp('......field map calculation')
    disp('.........saving quality map for masking')
    TE = (TEs(3)-TEs(1))*10^3; % effective echo time difference after phase complex fitting in seconds
    [~, FM_name,~] = fileparts(ph_V(1).fname) ;
    FM_romeo_file = sprintf('%s_romeo.nii',FM_name) ;
    if isunix
        status =  unix(sprintf('%s %s -m %s -o %s -t [%i,%i] -k nomask -g -q -B', para.romeo_command, ph_V(1).fname, mag_V(1).fname, FM_romeo_file, TE, TE)) ;
    elseif ispc
        status = system(sprintf('%s %s -m %s -o %s -t [%i,%i] -k nomask -g -q -B', para.romeo_command, ph_V(1).fname, mag_V(1).fname, FM_romeo_file, TE, TE)) ;
    end
    if status == 1
        error('ROMEO did not run properly - check your installation path')
    end
    
    
    
    %% field map rotation to scanner space
    % defining affine matrix in scanner space for data rotation to scanner
    % space with mantaining the same image origin (i.e. no translation)
    
    data_dim = size(ph_1tp.dat) ;
    Z = spm_imatrix(ph_1tp.mat) ;
    pixdim = Z(7:9);
    
    Maff_image = ph_1tp.mat ;
    O = Maff_image\[0 0 0 1]' ;
    O = O(1:3)' ;
    
    Maff_scanner(1,:) = [0 0 pixdim(3) -pixdim(3)*O(3)] ;
    Maff_scanner(2,:) = [pixdim(1) 0 0 -pixdim(1)*O(1)] ;
    Maff_scanner(3,:) = [0 pixdim(2) 0 -pixdim(2)*O(2)] ;
    Maff_scanner(4,:) = [0 0 0 1] ;
    
    
    FM = nifti('B0.nii') ;
    fm_data = FM.dat(:,:,:) ;
    fm_data(isnan(fm_data)) = 0 ;
    FM_V = spm_vol('B0.nii');
    spm_write_vol(FM_V, fm_data) ;
    clear FM
    img2scanner_mat = Maff_image\Maff_scanner ;
    FMrot = zeros(data_dim) ;
    
    for slice = 1 : data_dim(3)
        FMrot(:,:,slice) = spm_slice_vol(FM_V, img2scanner_mat*spm_matrix([0 0 slice]), data_dim(1:2), -7) ;
    end
    
    FMrot(isnan(FMrot)) = 0 ;
    FMrot_V = FM_V ;
    FMrot_V.mat = Maff_scanner ;
    FMrot_V.fname = 'B0_rot.nii';
    spm_write_vol(FMrot_V, FMrot)
    
    %% creating mask for QSM calculation
    
    disp('quality masking')
    qmask = nifti('quality.nii') ;
    mask_V = spm_vol('quality.nii') ;
    qmask = qmask.dat(:,:,:) ;
    qmask(isnan(qmask)) = 0 ;
    qmask(qmask>0.3) = 1 ;
    qmask(qmask<=0.3) = 0 ;
    qmask = imfill(qmask,6,'holes') ;
    qmask = smooth3(qmask, 'gaussian') ;
    qmask(qmask>0.6) = 1 ;
    qmask(qmask<=0.6) = 0 ;
    qmask = int16(qmask) ;
    mask_V.dt = [spm_type('int16') 0] ;
    mask_V.fname = 'mask.nii' ;
    spm_write_vol(mask_V, qmask) ;
    clear qmask
    qmask_rot = zeros(data_dim) ;
    for slice = 1 : data_dim(3)
        qmask_rot(:,:,slice) = spm_slice_vol(mask_V, img2scanner_mat*spm_matrix([0 0 slice]), data_dim(1:2), -7) ;
    end
    qmask_rot(isnan(qmask_rot)) = 0 ;
    mask_V.mat = Maff_scanner ;
    mask_V.fname = 'mask_rot.nii';
    spm_write_vol(mask_V, qmask_rot) ;
    
    
    
    %% SEPIA - background field removal and dipole inversion yielding final QSM
    
    disp('creating SEPIA header')
    B0 = para.B0 ;
    CF = B0*42.58*1e6;	% imaging frequency, in Hz (B0*gyromagnetic_ratio*1e6)
    delta_TE = 1;	    % echo spacing, in second - we have already combined data, in such situation set to 1
    B0_dir = [0;1;0];	% main magnetic field direction, it's always [0,1,0] because the images are resliced so that 2nd dimension is aligned with B0
    matrixSize = data_dim ;	    % image matrix size
    voxelSize = pixdim ;	% spatial resolution of the data, in mm
    header_fullfile = fullfile(output_fulldir, 'header_sepia.mat') ;
    save(header_fullfile, 'B0', 'B0_dir', 'CF', 'TE', 'delta_TE', 'matrixSize', 'voxelSize')
    
    % general SEPIA parameters
    sepia_addpath
    algorParam.general.isBET = 0 ;
    algorParam.general.isInvert = 1 ;
    algorParam.general.isGPU = 0 ;
    output_basename = fullfile(output_fulldir, 'sepia') ;
    
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
    algorParam.qsm.method = para.dipole_inv ;
    if strcmp(algorParam.qsm.method , 'ndi')
        
        algorParam.qsm.method = 'ndi' ;
        algorParam.qsm.tol = 1 ;
        algorParam.qsm.maxiter = 200 ;
        algorParam.qsm.stepSize = 1 ;
        
    elseif strcmp(algorParam.qsm.method , 'Star-QSM')
        
        algorParam.qsm.padsize = ones(1,3)*12 ;
        
    end
    
    
    disp('background field removal using PDF')
    BackgroundRemovalMacroIOWrapper(input,output_basename,input(2).name,algorParam);
    
    fprintf('dipole inversion using %s', algorParam.qsm.method)
    input(1).name = fullfile(output_fulldir, 'sepia_local-field.nii.gz') ;
    QSMMacroIOWrapper(input,output_basename,input(2).name,algorParam);
    
    
    disp('rotation of QSM back to the original image space')
    gunzip('sepia_QSM.nii.gz')
    QSM = nifti('sepia_QSM.nii') ;
    QSM = QSM.dat(:,:,:) ;
    
    scanner2img_mat = Maff_scanner\Maff_image ;
    QSMinvrot = zeros(data_dim) ;
    QSM_V = spm_vol('sepia_QSM.nii');
    for slice = 1 : data_dim(3)
        QSMinvrot(:,:,slice) = spm_slice_vol(QSM_V, scanner2img_mat*spm_matrix([0 0 slice]), data_dim(1:2), -7) ;
    end
    
    QSMinvrot_V = FM_V ;
    QSMinvrot_V.fname = 'sepia_QSM_invrot.nii';
    spm_write_vol(QSMinvrot_V, QSMinvrot)
    
    delete(sprintf('mag_TE%i.nii',size(TEs,2)));
    delete('sepia_mask-qsm.nii.gz')
    sprintf('finished after %s' , secs2hms(toc(tstart)))
end


