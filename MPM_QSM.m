% MPM QSM pipeline
% main steps:
% 1) complex-fit over echoes for pdw and t1w images
%    simple phase difference for mtw images
%    for odd and even echoes done separately
% 2) ROMEO phase unwrapping
% 3) masking based on ROMEO quality map
% 4) rotation to scanner space
% 5) PDF background field removal
% 6) star QSM for dpole inversion

% uses SEPIA toolbox
% Chan, K.-S., Marques, J.P., 2021. SEPIA—Susceptibility mapping pipeline tool for phase images. Neuroimage 227, 117611.
% script created by Barbara Dymerska
% last modifications 05/07/2021
% @ UCL FIL Physics


tstart = tic ;
%%%%% USER PARAMETERS %%%%%

% path to romeo phase uwnrapping followed by romeo command, i.e.
% (in linux) '/your_path/bin/romeo' or (in windows) 'D:\your_path\bin\romeo'
romeo_command = '~/Documents/MRI_software/ROMEO/romeo_linux_3.2.0/bin/romeo' ;

% for SEPIA header
B0 = 7;			    % magnetic field strength, in Tesla

% select dipole inversion method, either 'star' or 'ndi'
% 'ndi' may give more contrast but is less robust to noise
% 'Star-QSM' is very robust to noise and quick, may have less contrast than ndi
algorParam.qsm.method = 'Star-QSM' ;

% specify the angle of "slice rotation" (around 3rd dimention) in deg
alpha = 30 ;


% root directory to nifti files:
in_root_dir = '/media/barbara/hdd2/DATA/FIL/7T/20210623.M700198_FIL_analysis' ;
% the data will be saved in:
out_root_dir = '/media/barbara/hdd2/DATA/FIL/7T/20210623.M700198_FIL_analysis/SEPIA/MORSE_scan2';

% directories, parameters and files specific to given contrast:
for run = 1:3
    
        switch run
        case 1 %pdw
            mag_dir = 'pdw_mfc_3dflash_v1k_0025' ; % folder with magnitude niftis
            ph_dir = 'pdw_mfc_3dflash_v1k_0026' ; % folder with phase inftis
            TEs = [2.2 4.58 6.96 9.34 11.72 14.1] ; % echo time in ms
            output_dir = 'pdw_RR_25_26' ; % output directory for a specific submeasurement from MPM
            mag_file = 's2021-06-23_10-18-112654-00001-01728-6.nii' ; % magnitude reference nifti file for ROMEO unwrapping and masking
            
        case 2 % t1w
            mag_dir = 't1w_mfc_3dflash_v1k_0022' ;
            ph_dir = 't1w_mfc_3dflash_v1k_0023' ;
            TEs = [2.3 4.68 7.06 9.44 11.82 14.2] ; % echo time in ms
            output_dir = 't1w_RR_22_23' ;
            mag_file = 's2021-06-23_10-18-111631-00001-01728-6.nii' ; % magnitude reference nifti file for ROMEO unwrapping and masking
            
        case 3 % mtw
            mag_dir = 'mtw_mfc_3dflash_v1k_180deg_0031' ;
            ph_dir = 'mtw_mfc_3dflash_v1k_180deg_0032' ;
            TEs = [2.2 4.58 6.96 9.34] ; % echo time in ms
            output_dir = 'mtw_RR_31_32' ;
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
    
    
    clear ph_1tp.img mag_1tp
    
    % rescaling phase into [0,2*pi] phase range
    ph = 2*pi*single(ph - min(vector(ph)))/single(max(vector(ph))-min(vector(ph))) ;
    
    for read_dir = 1:2
        
        if size(TEs,2) < 6 % mtw has only 4 echoes so complex fit is not possible (minimum 3 echoes)
            disp('calculating phase difference')
            FM = angle(exp(1i*(ph(:,:,:,read_dir+2)-ph(:,:,:,read_dir)))) ;
        else
            disp('complex fitting phase')
            compl = single(mag).*exp(-1i*ph);
            [FM, dp1, relres, ~] = Fit_ppm_complex_TE(compl(:,:,:,read_dir:2:end),TEs(read_dir:2:end));
        end
        
        FM_both(:,:,:,read_dir) = FM ;
        
    end
    clear mag ph FM
    
    % saving odd and even echoes as one file for ROMEO unwrapping    
    FM_file = 'FM.nii' ;
    FM_both = make_nii(FM_both) ;
    FM_both.hdr.hist = ph_1tp.hdr.hist ;
    centre_and_save_nii(FM_both, FM_file, ph_1tp.hdr.dime.pixdim);
    clear FM_both
    
    % creating magnitude reference for ROMEO unwrapping
    mag_fullfile = fullfile(mag_fulldir, mag_file);
    mag_last = load_untouch_nii(mag_fullfile) ;
    mag_double = repmat(mag_last.img,[1 1 1 2]) ;
    make_nii(mag_double)
    mag_fullfile = sprintf('mag_TE%i.nii',size(TEs,2)) ;
    centre_and_save_nii(make_nii(mag_double), mag_fullfile, ph_1tp.hdr.dime.pixdim);
    clear mag_double
    
    disp('phase unwrapping with ROMEO + removing global mean value + field map calculation + saving quality map for masking')
    TE = (TEs(3)-TEs(1))*10^3; % effective echo time difference after phase complex fitting in seconds
    [~, FM_name,~] = fileparts(FM_file) ;
    FM_romeo_file = sprintf('%s_romeo.nii',FM_name) ;
    if isunix
        unix(sprintf('%s %s -m %s -o %s -t [%i,%i] -k nomask -g -q -B', romeo_command, FM_file, mag_fullfile, FM_romeo_file, TE, TE)) ;
    elseif ispc
        system(sprintf('%s %s -m %s -o %s -t [%i,%i] -k nomask -g -q -B', romeo_command, FM_file, mag_fullfile, FM_romeo_file, TE, TE)) ;
    end

    reslice_nii('B0.nii', 'B0_rot.nii',ph_1tp.hdr.dime.pixdim(2:4), 1, 0)
    FM = load_nii('B0_rot.nii') ;
    FM.img = changeImageSize(FM.img, circshift(ph_1tp.hdr.dime.dim(2:4),1)) ;
    FM.hdr.dime.dim(2:4) = circshift(ph_1tp.hdr.dime.dim(2:4),1) ;
    centre_and_save_nii(FM, 'B0_rot.nii', ph_1tp.hdr.dime.pixdim);
    
    
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
    qmask_file = fullfile(output_fulldir, 'mask.nii') ;
    save_nii(qmask, qmask_file);
    reslice_nii('mask.nii', 'mask_rot.nii', ph_1tp.hdr.dime.pixdim(2:4), 1 , 0)
    qmask_file = fullfile(output_fulldir, 'mask_rot.nii') ;
    qmask = load_nii(qmask_file) ;
    qmask.img = round(qmask.img) ;
    
    qmask.img = changeImageSize(qmask.img, circshift(ph_1tp.hdr.dime.dim(2:4),1)) ;
    qmask.hdr.dime.dim(2:4) = circshift(ph_1tp.hdr.dime.dim(2:4),1) ;
    centre_and_save_nii(qmask, qmask_file, ph_1tp.hdr.dime.pixdim);
    
    
    %% SEPIA - calculates QSM
    
    % create SEPIA header
    CF = B0*42.58*1e6;	% imaging frequency, in Hz (B0*gyromagnetic_ratio*1e6)
    delta_TE = 1;	    % echo spacing, in second - we have already combined data, in such situation set to 1
    B0_dir = [0;0;1];	% main magnetic field direction, it's always [0,0,1] because the images are resliced so that 3rd dimention is aligned with B0
    hdr = load_nii_hdr('FM_romeo_mean_rot.nii') ;
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
    algorParam.bfr.padSize = 40 ;
    
    
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
    
    sprintf('run %i finished after %s' ,run, secs2hms(toc))
    
    QSM = load_nii(fullfile(output_fulldir, 'QSM.nii.gz'));
    
    disp('rotation of QSM back to the original image space')
    QSM_invrot = flip(permute(QSM.img, [2 3 1]),1);
    
    M_aff(1,:) = [cos(deg2rad(alpha))    -sin(deg2rad(alpha))     0        0    ];
    M_aff(2,:) = [sin(deg2rad(alpha))    cos(deg2rad(alpha))      0        0    ];
    M_aff(3,:) = [0                      0                        1.0000   0    ];
    M_aff(4,:) = [0                      0                        0        1.0000];
    
    [QSM_invrot, ~] = affine(QSM_invrot, M_aff);
    QSM_invrot = changeImageSize(QSM_invrot, ph_1tp.hdr.dime.dim(2:4)) ;
    
    save_nii(make_nii(QSM_invrot), 'QSM_invrot.nii.gz')
    
    QSM_all(:,:,:,run) = QSM.img ;
    QSM_all_invrot(:,:,:,run) = QSM_invrot ;
    clear QSM QSM_invrot
    
end

QSM_all_mean = mean(QSM_all, 4) ;
QSM_pdw_t1w_mean = mean(QSM_all(:,:,:,1:2), 4) ;

QSM_all_invrot_mean = mean(QSM_all_invrot, 4) ;
QSM_pdw_t1w_invrot_mean = mean(QSM_all_invrot(:,:,:,1:2), 4) ;

centre_and_save_nii(make_nii(QSM_all_mean), fullfile(out_root_dir,'QSM_all_mean.nii'), ph_1tp.hdr.dime.pixdim);
centre_and_save_nii(make_nii(QSM_pdw_t1w_mean), fullfile(out_root_dir,'QSM_pdw_t1w_mean.nii'), ph_1tp.hdr.dime.pixdim);

centre_and_save_nii(make_nii(QSM_all_invrot_mean), fullfile(out_root_dir,'QSM_all_invrot_mean.nii'), ph_1tp.hdr.dime.pixdim);
centre_and_save_nii(make_nii(QSM_pdw_t1w_invrot_mean), fullfile(out_root_dir,'QSM_pdw_t1w_invrot_mean.nii'), ph_1tp.hdr.dime.pixdim);


sprintf('total processing finished after %s' , secs2hms(toc(tstart)))
