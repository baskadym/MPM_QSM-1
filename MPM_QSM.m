% MPM QSM pipeline
% it uses:
% 1) complex-fit over echoes for pdw and t1w images
%    simple phase difference for mtw images
%    for odd and even echoes separately
% 2) ROMEO phase unwrapping
% 3) SPM masking
% 4) PDF background field removal
% 5) non-linear dipole inversion
% uses SEPIA toolbox
% Chan, K.-S., Marques, J.P., 2021. SEPIA—Susceptibility mapping pipeline tool for phase images. Neuroimage 227, 117611.
% script created by Barbara Dymerska 22.06.2021
% @ UCL FIL Physics
%
% input variables:
% mag_dir - directory with magnitude niftis, one nifti per echoe
% ph_dir  - directory with phase niftis, one nifti per echoe
% mag4mask_file - selected magnitude for masking, usually last echo

%%%%% USER PARAMETERS %%%%%

mag4mask_file = '/media/barbara/hdd2/DATA/FIL/MORSE_Opt_phase/pdw_mfc_3dflash_v1k_RR_0054/s2021-05-24_14-26-153624-00001-01728-6.nii' ;
% masking parameters
thresh = 0.9 ;
fill_holes = 'yes' ;
conn = 8 ;

% for SEPIA header
%BKD: I do complex fit between odd and even echoes separately so the field map has effectively TE difference = TE(3)-TE(1)
B0 = 7;			    % magnetic field strength, in Tesla
B0_dir = [0;1;0];	% main magnetic field direction, [x,y,z]
CF = 7*42.58*1e6;	% imaging frequency, in Hz (B0*gyromagnetic_ratio*1e6)
delta_TE = 1;	    % echo spacing, in second - we have already combined data, in such situation set to 1
matrixSize = [364,426,288];	    % image matrix size
voxelSize = [0.6, 0.6, 0.6];	% spatial resolution of the data, in mm

% select dipole inversion method, either 'star' or 'ndi'
% 'ndi' may give more contrast but is less robust to noise
% 'star' is very robust to noise and quick, may have less contrast than ndi
algorParam.qsm.method = 'star' ;

for run = 1:3
    tic
    switch run
        
        case 1 %pdw
            mag_dir = '/media/barbara/hdd2/DATA/FIL/MORSE_Opt_phase/pdw_mfc_3dflash_v1k_RR_0054' ;
            ph_dir = '/media/barbara/hdd2/DATA/FIL/MORSE_Opt_phase/pdw_mfc_3dflash_v1k_RR_0055' ;
            TEs = [2.2 4.58 6.96 9.34 11.72 14.1] ; % echo time in ms
            output_dir = '/media/barbara/hdd2/DATA/FIL/MORSE_Opt_phase/SEPIA/pdw_RR_54_55/' ;
            
        case 2 % t1w
            mag_dir = '/media/barbara/hdd2/DATA/FIL/MORSE_Opt_phase/t1w_mfc_3dflash_v1k_RR_0056' ;
            ph_dir = '/media/barbara/hdd2/DATA/FIL/MORSE_Opt_phase/t1w_mfc_3dflash_v1k_RR_0057' ;
            TEs = [2.3 4.68 7.06 9.44 11.82 14.2] ; % echo time in ms
            output_dir = '/media/barbara/hdd2/DATA/FIL/MORSE_Opt_phase/SEPIA/t1w_RR_56_57/' ;
            
        case 3 % mtw
            mag_dir = '/media/barbara/hdd2/DATA/FIL/MORSE_Opt_phase/mtw_mfc_3dflash_v1k_180deg_RR_0058' ;
            ph_dir = '/media/barbara/hdd2/DATA/FIL/MORSE_Opt_phase/mtw_mfc_3dflash_v1k_180deg_RR_0059' ;
            TEs = [2.2 4.58 6.96 9.34] ; % echo time in ms
            output_dir = '/media/barbara/hdd2/DATA/FIL/MORSE_Opt_phase/SEPIA/mtw_RR_58_59/' ;
            
    end
    
    
    % path to romeo phase uwnrapping followed by romeo command, i.e.
    % (in linux) /executable/dir/romeo or (in windows) \\executable\dir\romeo
    romeo_command = '/media/barbara/hdd1/DATA/phase_unwrapping_comparison/ROMEO/romeo_linux_3.1.4/mnt/c/tmp/romeo_linux_3.1.4/bin/romeo' ;
    
    %%%%% END OF USER PARAMETERS %%%%%
    if ~exist(output_dir, 'dir')
        mkdir(output_dir)
    end
    cd(output_dir)
    
    TEs = TEs/10^3 ;
    
    for t = 1:size(TEs,2)
        file = dir(fullfile(ph_dir, sprintf('s20*-%i.nii', t)) );
        ph_1tp = load_untouch_nii(fullfile(ph_dir, file.name));
        ph(:,:,:,t) = ph_1tp.img ;
        
        if run ~= 3 % for mtw (run = 3) acquisition we cannot perform complex fit so we don't need magnitude data for each TE
            file = dir(fullfile(mag_dir, sprintf('s20*-%i.nii', t)) );
            mag_1tp = load_untouch_nii(fullfile(mag_dir, file.name)) ;
            mag(:,:,:,t) = mag_1tp.img ;
            
        end
        
    end
    
    
    clear ph_1tp.img mag_1tp
    
    % rescaling phase into [0,pi] phase range
    ph = 2*pi*single(ph)/4095 ;
    
    for read_dir = 1:2
        
        if run == 3 % mtw has only 4 echoes so complex fit is not possible (minimum 3 echoes)
            FM = angle(exp(1i*(ph(:,:,:,read_dir+2)-ph(:,:,:,read_dir)))) ;
        else % complex fit
            compl = single(mag).*exp(-1i*ph);
            [FM, ~, ~, ~] = Fit_ppm_complex_TE(compl(:,:,:,read_dir:2:end),TEs(read_dir:2:end));
            
        end
        
        if read_dir == 1
            flag = 'odd' ;
        else
            flag = 'even' ;
        end
        
        FM_file = fullfile(output_dir, sprintf('FM_%s.nii', flag)) ;
        centre_and_save_nii(make_nii(FM, ph_1tp.hdr.dime.pixdim(2:4)), FM_file , ph_1tp.hdr.dime.pixdim);
        
        % phase unwrapping with ROMEO, removing global mean as well
        [~, FM_name,~] = fileparts(FM_file) ;
        FM_romeo_file = fullfile(output_dir, sprintf('%s_romeo.nii',FM_name)) ;
        mag_file = dir(fullfile(mag_dir, sprintf('s20*-%i.nii', size(TEs,2))));
        unix(sprintf('%s -m %s -o %s -k nomask -g %s', romeo_command, fullfile(mag_file.folder, mag_file.name), FM_romeo_file, FM_file)) ;
    end
    
    
    % masking - mask from pdw acquisition is used for all 3 MPM acuisitions
    if run == 1
        [mask_file, ~] = mask_spm12_nofsl(mag4mask_file, thresh, fill_holes, conn) ;
    end
    
    % the echos are evenly spaces so I can simply average these after unwrapping and demeaning
    FM_odd = load_nii('FM_odd_romeo.nii') ;
    FM_even = load_nii('FM_even_romeo.nii') ;
    
    % averaging odd and even and scaling into the unit of Hz
    TE = (TEs(3)-TEs(1))/10^3; % effective echo time difference after phase complex fitting in seconds
    FM_mean = (FM_odd.img + FM_even.img)/(2*TE*2*pi) ;
    centre_and_save_nii(make_nii(FM_mean, ph_1tp.hdr.dime.pixdim(2:4)), 'FM_mean_romeo.nii' , ph_1tp.hdr.dime.pixdim);
    
    
    %% SEPIA - calculates QSM
    
    % create SEPIA header
    header_fullfile = fullfile(ph_dir, 'header_sepia.mat') ;
    save(header_fullfile, 'B0', 'B0_dir', 'CF', 'TE', 'delta_TE', 'matrixSize', 'voxelSize')
    
    % general SEPIA parameters
    sepia_addpath
    algorParam.general.isBET = 0 ;
    algorParam.general.isInvert = 1 ;
    algorParam.general.isGPU = 0 ;
    
    
    % inputs for background field removal
    
    input(1).name = fullfile(output_dir, 'FM_mean_romeo.nii') ;
    input(2).name = mask_file ;
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
        
    elseif strcmp(algorParam.qsm.method , 'star')
        
        algorParam.qsm.padsize = ones(1,3)*12 ;
        
    end
    
    
    output_basename = fullfile(output_dir, sprintf('sepia_%s_%s', algorParam.bfr.method, algorParam.qsm.method)) ;
    sprintf('run %i preprocessing finished after %s' ,run, secs2hms(toc))
    
    % background field removal
    BackgroundRemovalMacroIOWrapper(input,output_basename,mask_file,algorParam);
    
    % dipole inversion
    input(1).name = fullfile(output_dir, sprintf('sepia_%s_%s_local-field.nii.gz', algorParam.bfr.method, algorParam.qsm.method)) ;
    QSMMacroIOWrapper(input,output_basename,mask_file,algorParam);
    
    sprintf('run %i finished after %s' ,run, secs2hms(toc))
    
end
exit