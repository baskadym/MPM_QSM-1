function [mask_file, mask] = mask_spm12(mag_file, thresh, fill_holes, conn, mask_mag)
% Barbara Dymerska
% brain extraction script using spm12 segment tool
% it sums three probability maps from segment tool, thresholds the resulting map to obtain binary mask and fills residual holes
% Inputs:
% mag_file      - file from which to estimate the mask
% thresh        - between 0 and 1, smaller value bigger mask, threshold of the spm probability maps to create a binary map
% fill_holes    - true or false, if true it will fill small holes within the mask
% conn          - used in imfill.m, connectivity of neighbouring voxels, choose one of the following: 4, 6, 8, 18, 26
% mask_mag      - true or false, if true it will also output masked mag_file

tic
%% calling SPM segment tool
spm('defaults', 'FMRI');
matlabbatch{1}.spm.spatial.preproc.channel.vols = {sprintf('%s,1',mag_file)};
matlabbatch{1}.spm.spatial.preproc.channel.biasreg = 0.001;
matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 60;
matlabbatch{1}.spm.spatial.preproc.channel.write = [1 1];
matlabbatch{1}.spm.spatial.preproc.tissue(1).tpm = {'/home/bdymerska/Documents/github/SPM/spm12/tpm/TPM.nii,1'};
matlabbatch{1}.spm.spatial.preproc.tissue(1).ngaus = 1;
matlabbatch{1}.spm.spatial.preproc.tissue(1).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(1).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(2).tpm = {'/home/bdymerska/Documents/github/SPM/spm12/tpm/TPM.nii,2'};
matlabbatch{1}.spm.spatial.preproc.tissue(2).ngaus = 1;
matlabbatch{1}.spm.spatial.preproc.tissue(2).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(2).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(3).tpm = {'/home/bdymerska/Documents/github/SPM/spm12/tpm/TPM.nii,3'};
matlabbatch{1}.spm.spatial.preproc.tissue(3).ngaus = 2;
matlabbatch{1}.spm.spatial.preproc.tissue(3).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(3).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(4).tpm = {'/home/bdymerska/Documents/github/SPM/spm12/tpm/TPM.nii,4'};
matlabbatch{1}.spm.spatial.preproc.tissue(4).ngaus = 3;
matlabbatch{1}.spm.spatial.preproc.tissue(4).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(4).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(5).tpm = {'/home/bdymerska/Documents/github/SPM/spm12/tpm/TPM.nii,5'};
matlabbatch{1}.spm.spatial.preproc.tissue(5).ngaus = 4;
matlabbatch{1}.spm.spatial.preproc.tissue(5).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(5).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(6).tpm = {'/home/bdymerska/Documents/github/SPM/spm12/tpm/TPM.nii,6'};
matlabbatch{1}.spm.spatial.preproc.tissue(6).ngaus = 2;
matlabbatch{1}.spm.spatial.preproc.tissue(6).native = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(6).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.warp.mrf = 1;
matlabbatch{1}.spm.spatial.preproc.warp.cleanup = 1;
matlabbatch{1}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
matlabbatch{1}.spm.spatial.preproc.warp.affreg = 'mni';
matlabbatch{1}.spm.spatial.preproc.warp.fwhm = 0;
matlabbatch{1}.spm.spatial.preproc.warp.samp = 3;
matlabbatch{1}.spm.spatial.preproc.warp.write = [0 0];

spm_jobman('serial',matlabbatch) ;

%% summing up probability maps for GM, WM and CSF
[path,mag_name,~] = fileparts(mag_file) ;
out_dir = fullfile(path, 'spm_mask');
mkdir(out_dir)
mask_file = fullfile(out_dir, 'mask.nii') ;

c1 = nifti(fullfile(path, sprintf('c1%s.nii', mag_name))) ;
c2 = nifti(fullfile(path, sprintf('c2%s.nii', mag_name))) ;
c3 = nifti(fullfile(path, sprintf('c3%s.nii', mag_name))) ;

mask(:,:,:) = c1.dat(:,:,:) + c2.dat(:,:,:) + c3.dat(:,:,:) ;
mask(mask>thresh) = 1 ;
mask(mask<=thresh) = 0 ;

if fill_holes
    mask = imfill(mask,conn,'holes') ;
else
    warning('not filling holes')
end

createNifti(mask, mask_file, c1.mat) ;

if mask_mag
    mag_masked_obj = nifti(mag_file) ;
    mag_masked = mag_masked_obj.dat(:,:,:).*mask(:,:,:) ;
    mag_masked_obj.dat.fname = fullfile(out_dir, sprintf('%s_masked.nii',mag_name));
    create(mag_masked_obj)
    mag_masked_obj.dat(:,:,:) = mag_masked ;
end

for j = 1:5
movefile(fullfile(path,sprintf('c%i%s.nii',j,mag_name)), out_dir)
end
movefile(fullfile(path,sprintf('m%s.nii',mag_name)), out_dir)
movefile(fullfile(path,sprintf('%s_seg8.mat',mag_name)), out_dir)
movefile(fullfile(path,sprintf('BiasField_%s.nii',mag_name)), out_dir)

disp(['spm masking took: ' secs2hms(toc)]);

