% this is the script used for voxel based analysis

% dataset name
%dataset = '/indirect/staff/cyrilpernet/multispectral_segmentation/ds003653'; 
dataset = '/indirect/staff/cyrilpernet/multispectral_segmentation/nrudataset'; 

% set the SPM parameters
global defaults
defaults = spm_get_defaults;
defaults.modality= 'PET'; % makes no difference. that's just what SPM does for VBM
defaults.stats.rft.nonstat = 1; % matters for RFT 
spm_jobman('initcfg')

% make and run the batch

% destination is NOT in the shareable folder, just to big for github
% the maps are however on NeuroVault
destination = [dataset filesep 'derivatives' filesep 'RepANOVA'];
mkdir(destination)
matlabbatch{1}.spm.stats.factorial_design.dir = {destination};

matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(1).name     = 'subject';
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(1).dept     = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(1).variance = 1;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(1).gmsca    = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(1).ancova   = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(2).name     = 'input';
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(2).dept     = 1;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(2).variance = 1;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(2).gmsca    = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(2).ancova   = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(3).name     = 'model';
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(3).dept     = 1;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(3).variance = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(3).gmsca    = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(3).ancova   = 0;

% get all the data ; from running the 'segment_images.m we know the
% subfolders where to find the data
cond1 = dir(fullfile(dataset,'derivatives/T1_nG1/wc1*.nii'));
cond2 = dir(fullfile(dataset,'derivatives/T1_nG2/wc1*.nii'));
cond3 = dir(fullfile(dataset,'derivatives/T12_nG1/wc1*.nii'));
cond4 = dir(fullfile(dataset,'derivatives/T12_nG2/wc1*.nii'));
if size(cond1,1) ~= size(cond2,1) || ...
        size(cond2,1) ~= size(cond3,1) || ...
        size(cond3,1) ~= size(cond4,1)
    error('conditions (i.e. segmentation types) do not have the same nu,ber of subjects')
end

for s=1:size(cond1,1)
    scans{1} = [fullfile(cond1(s).folder,cond1(s).name) ',1'];
    scans{2} = [fullfile(cond2(s).folder,cond2(s).name) ',1'];
    scans{3} = [fullfile(cond3(s).folder,cond3(s).name) ',1'];
    scans{4} = [fullfile(cond4(s).folder,cond4(s).name) ',1'];
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fsuball.fsubject(s).scans = scans';
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fsuball.fsubject(s).conds = ...
        [1 1
        1 2
        2 1
        2 2];
end

matlabbatch{1}.spm.stats.factorial_design.des.fblock.maininters{1}.fmain.fnum = 1;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.maininters{2}.fmain.fnum = 2;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.maininters{3}.inter.fnums = [1 2];

matlabbatch{1}.spm.stats.factorial_design.cov       = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});

matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none     = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im             = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em             = {''};

matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit         = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm        = 1;

matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

% run
out=spm_jobman('run',matlabbatch);