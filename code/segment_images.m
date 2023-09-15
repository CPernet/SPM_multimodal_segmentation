function out = segment_images(datadir,outdir,options)
% perform the SPM12 segmentation, get tissue volume information, get voxel
% distributions of normalized images, compute the DARTEL template and
% return the jobs
%
% INPUTS datadir directory where all the data are
%        outdir where to save the volumes and distrib mat files
%        options structure with 'modality', 'NGaussian','Modulate' options
%        debug internal checks
%
% OUPTPUT out is a cell array of the segmentation + dartel jobs from SPM job manager
%
% Cyril Pernet - CCBS, University of Edinburgh / NRU, Rigshospitalet (2022) - CCBY
% Modifed by Marc Cummings - 2023

%% basic info
%-----------------------------------------------------------------------
cd(datadir)
spmroot = fileparts(which('spm'));
spm('defaults', 'FMRI'); spm_jobman('initcfg')
BIDS = spm_BIDS(pwd);

%% create batch per subjects
%-----------------------------------------------------------------------
index = 1;
prefileT1 = [spm_BIDS(BIDS,'data','type','T1w')];
prefileT2 = [spm_BIDS(BIDS,'data','type','T2w')];
% remove compressed file if uncompressed files exist
if (all(contains(prefileT1,'.nii.gz')) && all(contains(prefileT2,'.nii.gz')))
    prefileT1((find(endsWith(prefileT1,'.nii') == 1) + 1)) = [];
    prefileT2((find(endsWith(prefileT2,'.nii') == 1) + 1)) = [];
end
% create file Map over all the files of the T1 and T2
fileMap = struct('ID',{},'path',{},'T1Path',{},'T2Path',{});
for i = 1:length(prefileT1)
    fileMap(i).ID = BIDS.subjects(i).name;
    fileMap(i).path = append(BIDS.subjects(i).path, filesep, 'anat', filesep);
    fileMap(i).T1Path = prefileT1(i);
    fileMap(i).T2Path = prefileT2(i);
end

for mapIndex = length(fileMap):-1:1
    subjectFiles           = fileMap(mapIndex);
    subjectFiles.T1Path{1} = decompress_gzip(subjectFiles.T1Path{1});

    % cleanup
    [filepath,name,ext] = fileparts(subjectFiles.T1Path{1});
    T1name              = append(name, ext);
    c    = dir([filepath filesep 'c*' T1name]);  for d=1:length(c); delete(fullfile(c(d).folder,c(d).name)); end
    rc   = dir([filepath filesep 'rc*' T1name]); for d=1:length(rc); delete(fullfile(rc(d).folder,rc(d).name)); end
    wc   = dir([filepath filesep 'wc*' T1name]); for d=1:length(wc); delete(fullfile(wc(d).folder,wc(d).name)); end
    seg  = dir(append(filepath, filesep, '*seg8.mat')); if ~isempty(seg); delete(fullfile(seg.folder,seg.name)); end
    u_rc = dir(append(filepath, filesep, 'u_rc*')); if ~isempty(u_rc); delete(fullfile(u_rc.folder,u_rc.name)); end

    % channel T1 or T1 and coregistered T2
    if strcmpi(options.modality,'T1')
        matlabbatch{1}.spm.spatial.preproc.channel(1).vols(1) = subjectFiles.T1Path;
        matlabbatch{1}.spm.spatial.preproc.channel(1).biasreg = 0.001;
        matlabbatch{1}.spm.spatial.preproc.channel(1).biasfwhm = 60;
        matlabbatch{1}.spm.spatial.preproc.channel(1).write = [0 0];
        batch_index = 1;

    elseif strcmpi(options.modality,'T12')
        matlabbatch{1}.spm.spatial.coreg.estwrite.ref(1) = subjectFiles.T1Path;
        subjectFiles.T2Path{1} = decompress_gzip(subjectFiles.T2Path{1});
        matlabbatch{1}.spm.spatial.coreg.estwrite.source = subjectFiles.T2Path;
        matlabbatch{1}.spm.spatial.coreg.estwrite.other = {''};
        matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
        matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
        matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
        matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
        matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 4;
        matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
        matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
        matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';

        matlabbatch{2}.spm.spatial.preproc.channel(1).vols(1) = subjectFiles.T1Path;
        matlabbatch{2}.spm.spatial.preproc.channel(1).biasreg = 0.001;
        matlabbatch{2}.spm.spatial.preproc.channel(1).biasfwhm = 60;
        matlabbatch{2}.spm.spatial.preproc.channel(1).write = [0 0];
        matlabbatch{2}.spm.spatial.preproc.channel(2).vols(1) = cfg_dep('Coregister: Estimate & Reslice: Resliced Images', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','rfiles'));
        matlabbatch{2}.spm.spatial.preproc.channel(2).biasreg = 0.001;
        matlabbatch{2}.spm.spatial.preproc.channel(2).biasfwhm = 60;
        matlabbatch{2}.spm.spatial.preproc.channel(2).write = [0 0];
        batch_index = 2;
    end

    % specify nb of Gaussians and Modulations
    if strcmp(options.Modulate,'Yes')
        warped = [1 1];
    else
        warped = [1 0];
    end
    matlabbatch{batch_index}.spm.spatial.preproc.tissue(1).tpm = {[spmroot filesep 'tpm' filesep 'TPM.nii,1']};
    matlabbatch{batch_index}.spm.spatial.preproc.tissue(1).ngaus = options.NGaussian;
    matlabbatch{batch_index}.spm.spatial.preproc.tissue(1).native = [0 1];
    matlabbatch{batch_index}.spm.spatial.preproc.tissue(1).warped = warped;
    matlabbatch{batch_index}.spm.spatial.preproc.tissue(2).tpm = {[spmroot filesep 'tpm' filesep 'TPM.nii,2']};
    matlabbatch{batch_index}.spm.spatial.preproc.tissue(2).ngaus = options.NGaussian;
    matlabbatch{batch_index}.spm.spatial.preproc.tissue(2).native = [0 1];
    matlabbatch{batch_index}.spm.spatial.preproc.tissue(2).warped = warped;
    matlabbatch{batch_index}.spm.spatial.preproc.tissue(3).tpm = {[spmroot filesep 'tpm' filesep 'TPM.nii,3']};
    matlabbatch{batch_index}.spm.spatial.preproc.tissue(3).ngaus = options.NGaussian;
    matlabbatch{batch_index}.spm.spatial.preproc.tissue(3).native = [0 1];
    matlabbatch{batch_index}.spm.spatial.preproc.tissue(3).warped = warped;
    matlabbatch{batch_index}.spm.spatial.preproc.tissue(4).tpm = {[spmroot filesep 'tpm' filesep 'TPM.nii,4']};
    matlabbatch{batch_index}.spm.spatial.preproc.tissue(4).ngaus = 3;
    matlabbatch{batch_index}.spm.spatial.preproc.tissue(4).native = [0 1];
    matlabbatch{batch_index}.spm.spatial.preproc.tissue(4).warped = [1 0];
    matlabbatch{batch_index}.spm.spatial.preproc.tissue(5).tpm = {[spmroot filesep 'tpm' filesep 'TPM.nii,5']};
    matlabbatch{batch_index}.spm.spatial.preproc.tissue(5).ngaus = 4;
    matlabbatch{batch_index}.spm.spatial.preproc.tissue(5).native = [0 1];
    matlabbatch{batch_index}.spm.spatial.preproc.tissue(5).warped = [0 0];
    matlabbatch{batch_index}.spm.spatial.preproc.tissue(6).tpm = {[spmroot filesep 'tpm' filesep 'TPM.nii,6']};
    matlabbatch{batch_index}.spm.spatial.preproc.tissue(6).ngaus = 2;
    matlabbatch{batch_index}.spm.spatial.preproc.tissue(6).native = [0 1];
    matlabbatch{batch_index}.spm.spatial.preproc.tissue(6).warped = [0 0];
    matlabbatch{batch_index}.spm.spatial.preproc.warp.mrf = 1;
    matlabbatch{batch_index}.spm.spatial.preproc.warp.cleanup = 1;
    matlabbatch{batch_index}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
    matlabbatch{batch_index}.spm.spatial.preproc.warp.affreg = 'mni';
    matlabbatch{batch_index}.spm.spatial.preproc.warp.fwhm = 0;
    matlabbatch{batch_index}.spm.spatial.preproc.warp.samp = 3;
    matlabbatch{batch_index}.spm.spatial.preproc.warp.write = [0 0];

    % compute volumes
    matlabbatch{batch_index+1}.spm.util.tvol.matfiles(1) = cfg_dep('Segment: Seg Params', substruct('.','val', '{}',{batch_index}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','param', '()',{':'}));
    matlabbatch{batch_index+1}.spm.util.tvol.tmax = 6;
    matlabbatch{batch_index+1}.spm.util.tvol.mask = {[spmroot filesep 'tpm' filesep 'mask_ICV.nii,1']};
    matlabbatch{batch_index+1}.spm.util.tvol.outf = T1name;

    % save into a 'batch' variable
    batch{index} = matlabbatch;
    index = index+1;
    clear matlabbatch
end

%% run the jobs in parallel
%-----------------------------------------------------------------------
Ncores = 12; % set to use 12 processors in parallel
p      = gcp('nocreate');
if isempty(p)
    c            = parcluster;
    c.NumWorkers = Ncores;
    saveProfile(c);
    try
        parpool(Ncores-1);
    catch flasestart %#ok<NASGU>
        delete(c.Jobs)
        parcluster('local')
        parpool(Ncores-1);
    end
end

N = length(batch);
spm_jobman('initcfg')
parfor subject=1:N
    try
        % spm_jobman('initcfg') - depends server settings some need that
        out{subject}=spm_jobman('run',batch{subject});
    catch
        out{subject} = sprintf('errorbatch %g',subject);
    end
end
% delete(gcp('nocreate')) % close parallel pool

%% for each images get tissue distributions and metrics
%-----------------------------------------------------------------------
for p=4:-1:1
    P(p,:) = [spmroot filesep 'tpm' filesep 'TPM.nii,' num2str(p)]; % template 1 to 4
end
V           = spm_vol(P); clear P;
M           = single(spm_read_vols(V));
M           = mean(M,4)>0.01;            % mask image (GM,WM,CSF,Meninges)
[x,y,z]     = ind2sub(V(1).dim,find(M)); % location of these voxels

% Load vessels and nuclei templates
vessels_template = [fileparts(which('multispectral_segmentation_analysis.m')), ...
    filesep, 'Atlases', filesep, 'Vessels', filesep, ...
    'rvesselRadius.nii'];
nuclei_template   = [fileparts(which('multispectral_segmentation_analysis.m')), ...
    filesep, 'Atlases', filesep, 'Nuclei', filesep, ...
    'rprob_atlas_bilateral.nii'];

% Load template volumes and create masks
vessels_vol   = spm_vol(vessels_template);
vessels_mask  = spm_read_vols(vessels_vol) > 0.5; % diameters of vessels > 0.5 mm
nuclei_vol    = spm_vol(nuclei_template);
nuclei_mask   = spm_read_vols(nuclei_vol) > 0;

% Find voxel coordinates for the masks and store in variables
[vessels_x, vessels_y, vessels_z] = ind2sub(vessels_vol.dim, find(vessels_mask));
distrib_vessels                   = NaN(length(vessels_x),N,3);

for v = size(nuclei_vol,1):-1:1
    [nuclei_x{v}, nuclei_y{v}, nuclei_z{v}] = ind2sub(nuclei_vol(v).dim, find(squeeze(nuclei_mask(:,:,:,v))));
    distrib_nuclei{v}                       = NaN(length(nuclei_x{v}),N,3);
end

% compute
for subject=N:-1:1

    % get the volume information in ml
    if isfield(out{subject}{2},'vol1')
        volumes(subject,:) = [out{subject}{2}.vol1, out{subject}{2}.vol2,...
            out{subject}{2}.vol3,out{subject}{2}.vol4,...
            out{subject}{2}.vol5,out{subject}{2}.vol6];
    else
        volumes(subject,:) = [out{subject}{3}.vol1, out{subject}{3}.vol2,...
            out{subject}{3}.vol3,out{subject}{3}.vol4,...
            out{subject}{3}.vol5,out{subject}{3}.vol6];
    end

    % create array to hold tissue voxels
    tmpGM       = spm_read_vols(spm_vol(fullfile(spmroot,['tpm' filesep 'TPM.nii,1'])));
    tmp_tissues = NaN(size(tmpGM,1), size(tmpGM,2), size(tmpGM,3), 3);

    % get the in mask voxel distributions
    for tissue_class = 1:3
        tmp = dir([fileMap(subject).path filesep 'wc' num2str(tissue_class) '*.nii']);
        % copy file for voxelwise analysis
        if  ~exist([outdir filesep options.modality '_nG' num2str(options.NGaussian)],'dir')
            mkdir([outdir filesep options.modality '_nG' num2str(options.NGaussian)])
        end
        copyfile(fullfile(tmp.folder, tmp.name),...
            fullfile([outdir filesep options.modality '_nG' num2str(options.NGaussian)], tmp.name))
        distrib(:,subject,tissue_class)               = spm_get_data(fullfile(tmp.folder, tmp.name),[x,y,z]');
        tmp_tissues(:,:,:,tissue_class)               = spm_read_vols(spm_vol(fullfile(tmp.folder, tmp.name)));
        distrib_vessels(:,subject,tissue_class)       = spm_get_data(fullfile(tmp.folder, tmp.name), [vessels_x, vessels_y, vessels_z]');
        for v = size(nuclei_vol,1):-1:1
            distrib_nuclei{v}(:,subject,tissue_class) = spm_get_data(fullfile(tmp.folder, tmp.name), [nuclei_x{v}, nuclei_y{v}, nuclei_z{v}]');
        end

        % calculate entropy for tissue
        entropy(subject,tissue_class) = image_entropy(fullfile(tmp.folder, tmp.name));
    end

    % calculate the dunn index for tissues
    dunnIndexes(subject,:) = modified_DunnIndex(squeeze(tmp_tissues(:,:,:,1)),....
        squeeze(tmp_tissues(:,:,:,2)),squeeze(tmp_tissues(:,:,:,3)));

end
% save information
temp_name = ['volumes' options.modality '_nG' num2str(options.NGaussian)];
save(fullfile(outdir, temp_name), 'volumes', '-v7.3')
temp_name = ['distrib' options.modality '_nG' num2str(options.NGaussian)];
save(fullfile(outdir, temp_name), 'distrib', '-v7.3')
temp_name = ['dunnIndex' options.modality '_nG' num2str(options.NGaussian)];
save(fullfile(outdir, temp_name), 'dunnIndexes', '-v7.3')
temp_name = ['entropy' options.modality '_nG' num2str(options.NGaussian)];
save(fullfile(outdir, temp_name), 'entropy', '-v7.3')
temp_name = ['distrib_nuclei' options.modality '_nG' num2str(options.NGaussian)];
save(fullfile(outdir, temp_name), 'distrib_nuclei', '-v7.3')
temp_name = ['distrib_vessels' options.modality '_nG' num2str(options.NGaussian)];
save(fullfile(outdir, temp_name), 'distrib_vessels', '-v7.3')
clear distrib volumes dunnIndexes entropy tmpGM tmp_tissues distrib_nuclei distrib_vessels

%% generate the DARTEL template
%-----------------------------------------------------------------------

for subject=N:-1:1
    for tissue_class = 6:-1:1
        tmp = dir([fileMap(subject).path filesep 'rc' num2str(tissue_class) '*.nii']);
        if tissue_class == 1
            c1{subject} = [tmp.folder filesep tmp.name];
        elseif tissue_class == 2
            c2{subject} = [tmp.folder filesep tmp.name];
        elseif tissue_class == 3
            c3{subject} = [tmp.folder filesep tmp.name];
        elseif tissue_class == 4
            c4{subject} = [tmp.folder filesep tmp.name];
        elseif tissue_class == 5
            c5{subject} = [tmp.folder filesep tmp.name];
        elseif tissue_class == 6
            c6{subject} = [tmp.folder filesep tmp.name];
        end
    end
end

cd(datadir)
matlabbatch{1}.spm.tools.dartel.warp.images = {c1',c2',c3',c4',c5',c6'};
temp_name = ['template' options.modality '_nG' num2str(options.NGaussian)];
matlabbatch{1}.spm.tools.dartel.warp.settings.template = temp_name;
matlabbatch{1}.spm.tools.dartel.warp.settings.rform = 0;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(1).its = 3;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(1).rparam = [4 2 1e-06];
matlabbatch{1}.spm.tools.dartel.warp.settings.param(1).K = 0;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(1).slam = 16;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(2).its = 3;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(2).rparam = [2 1 1e-06];
matlabbatch{1}.spm.tools.dartel.warp.settings.param(2).K = 0;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(2).slam = 8;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(3).its = 3;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(3).rparam = [1 0.5 1e-06];
matlabbatch{1}.spm.tools.dartel.warp.settings.param(3).K = 1;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(3).slam = 4;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(4).its = 3;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(4).rparam = [0.5 0.25 1e-06];
matlabbatch{1}.spm.tools.dartel.warp.settings.param(4).K = 2;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(4).slam = 2;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(5).its = 3;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(5).rparam = [0.25 0.125 1e-06];
matlabbatch{1}.spm.tools.dartel.warp.settings.param(5).K = 4;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(5).slam = 1;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(6).its = 3;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(6).rparam = [0.25 0.125 1e-06];
matlabbatch{1}.spm.tools.dartel.warp.settings.param(6).K = 6;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(6).slam = 0.5;
matlabbatch{1}.spm.tools.dartel.warp.settings.optim.lmreg = 0.01;
matlabbatch{1}.spm.tools.dartel.warp.settings.optim.cyc = 3;
matlabbatch{1}.spm.tools.dartel.warp.settings.optim.its = 3;
out{length(out)+1} = spm_jobman('run',matlabbatch);


