function out = segment_images(datadir,outdir,options)
% perform the SPM12 segmentation, get tissue volume information, get voxel
% distributions of normalized images, compute the DARTEL template and
% return the jobs
%
% INPUTS datadir directory where all the data are 
%                 this is not BIDS so the mapping is done in the code (i.e. either map data as below or recode, sorry)
%                 this is one folder per subject with subdirectory 'newT1' and 'T2'
%                     sub-001/newT1/somename.nii
%                     sub-001/T2/somename.nii
%                     sub-002/newT1/somename.nii
%                     sub-002/T2/somename.nii
%        outdir where to save the volumes and distrib mat files
%        options structure with 'modality', 'NGaussian','Modulate' options
%
% OUPTPUT out is a cell array of the segmentation + dartel jobs from SPM job manager
%
% Cyril Pernet - CCBS, University of Edinburgh / NRU, Rigshospitalet (2022) - CCBY

%% basic info
%-----------------------------------------------------------------------
cd(datadir)
folders = dir; % maps the subjects' folders
spmroot = fileparts(which('spm'));
spm('defaults', 'FMRI'); spm_jobman('initcfg')

%% create batch per subjects
%-----------------------------------------------------------------------
index = 1;
for subject = 3:length(folders)
    if folders(subject).isdir
        
       cd(fullfile(datadir, folders(subject).name));
       name = folders(subject).name;
       folders = dir; % maps the sessions' folders
       for session = 3:length(folders)
            if folders(session).isdir
                name = name + "_" + folders(session).name;
                cd(folders(session).name + "/anat");
                decompress_gzip();
                % channel T1 or T1 and coregistered T2
                if strcmpi(options.modality,'T1')
                    T1name = dir('*_T1w.nii');
                    c   = dir(['c*' T1name.name]);  for d=1:length(c); delete(fullfile(c(d).folder,c(d).name)); end
                    rc  = dir(['rc*' T1name.name]); for d=1:length(rc); delete(fullfile(rc(d).folder,rc(d).name)); end
                    wc  = dir(['wc*' T1name.name]); for d=1:length(wc); delete(fullfile(wc(d).folder,wc(d).name)); end
                    seg = dir('*seg8.mat'); if ~isempty(seg); delete(fullfile(seg.folder,seg.name)); end
                    u_rc = dir('u_rc*'); if ~isempty(u_rc); delete(fullfile(u_rc.folder,u_rc.name)); end
                    T1name = [pwd filesep T1name.name];
                    matlabbatch{1}.spm.spatial.preproc.channel(1).vols(1) = {T1name};
                    matlabbatch{1}.spm.spatial.preproc.channel(1).biasreg = 0.001;
                    matlabbatch{1}.spm.spatial.preproc.channel(1).biasfwhm = 60;
                    matlabbatch{1}.spm.spatial.preproc.channel(1).write = [0 0];
                    batch_index = 1;
                elseif strcmpi(options.modality,'T12')
                    T1name = dir('*_T1w.nii');
                    T1name = [pwd filesep T1name.name];
                    [pth, T1name, ext] = fileparts(T1name);
                    T1name = [T1name ext];
                    
                    cd(['..' filesep 'T2w']);
                    T2name = dir('*_T2w.nii');
                    N = length(T2name);
                    for n=1:N
                        try
                            if strcmp(T2name(n).name,T1name)
                                T2name(n) = [];
                            end
                        end
                    end
                    c = dir(['c*' T1name]); for d=1:length(c); delete(fullfile(c(d).folder,c(d).name)); end
                    rc = dir(['rc*' T1name]); for d=1:length(rc); delete(fullfile(rc(d).folder,rc(d).name)); end
                    wc = dir(['wc*' T1name]); for d=1:length(wc); delete(fullfile(wc(d).folder,wc(d).name)); end
                    seg= dir('*seg8.mat'); if ~isempty(seg); delete(fullfile(seg.folder,seg.name)); end
                    u_rc = dir('u_rc*'); if ~isempty(u_rc); delete(fullfile(u_rc.folder,u_rc.name)); end
                    T1name = [pwd filesep T1name];
                    T2name = [pwd filesep T2name(1).name];
                    
                    matlabbatch{1}.spm.spatial.coreg.estwrite.ref(1) = {T1name};
                    matlabbatch{1}.spm.spatial.coreg.estwrite.source = {T2name};
                    matlabbatch{1}.spm.spatial.coreg.estwrite.other = {''};
                    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
                    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
                    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
                    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
                    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 4;
                    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
                    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
                    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';
                    
                    matlabbatch{2}.spm.spatial.preproc.channel(1).vols(1) = {T1name};
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
                matlabbatch{batch_index+1}.spm.util.tvol.tmax = 3;
                matlabbatch{batch_index+1}.spm.util.tvol.mask = {[spmroot filesep 'tpm' filesep 'mask_ICV.nii,1']};
                matlabbatch{batch_index+1}.spm.util.tvol.outf = [fileparts(T1name)];
                
                % save into a 'batch' variable
                batch{index} = matlabbatch;
                index = index+1;
                clear matlabbatch
            end
       end
    end
end

%% run the jobs in parallel
%-----------------------------------------------------------------------
disp("run the jobs in parallel")
cd(datadir)
N = length(batch);
for subject=1:N
    try
        % spm_jobman('initcfg')
        out{subject}=spm_jobman('run',batch{subject});
    catch
        out{subject} = sprintf('errorbatch %g',subject);
    end
end


%% for each images get the metadata
%-----------------------------------------------------------------------
for p=1:4
    P(p,:) = [spmroot filesep 'tpm' filesep 'TPM.nii,' num2str(p)]; % template 1 to 4
end
V          = spm_vol(P);
M          = single(spm_read_vols(V));
M          = mean(M,4)>0.01;            % mask image (GM,WM,CSF,Meninges)
nvoxels    = sum(M(:));                 % how many in mask voxels
[x,y,z]    = ind2sub(V(1).dim,find(M)); % location of these voxels
distrib    = NaN(nvoxels,N,3);          % matrix of all voxels by N subjects by 3 tissue classes
volumes    = NaN(N,3);                  % matrix of N subjects by 3 tissue classes

for subject=3:(N+2)
    
    % get the volume information in ml
    cd(datadir)
    cd(folders(subject).name + "/anat");
    
    try
        f = dir('*seg8.mat'); results = load(f.name);
        disp(results);
        volumes(subject-2,:) = results.volumes.litres*1000;
        
        % get the in mask voxel distributions
        for tissue_class = 1:3
            tmp = dir(['wc' num2str(tissue_class) '*.nii']);
            distrib(:,subject-2,tissue_class) = spm_get_data(tmp.name,[x,y,z]');
        end
    end
end
temp_name = ['volumes' options.modality '_nG' num2str(options.NGaussian)];
save(outdir+filesep+temp_name,'volumes','-v7.3')
temp_name = ['distrib' options.modality '_nG' num2str(options.NGaussian)];
save(outdir+filesep+temp_name,'distrib','-v7.3')
clear distrib volumes

%% generate the DARTEL template
%-----------------------------------------------------------------------

% for subject=1:N
%     c1{subject} = cell2mat(out{subject}{1}.tiss(1).rc);
%     c2{subject} = cell2mat(out{subject}{1}.tiss(2).rc);
%     c3{subject} = cell2mat(out{subject}{1}.tiss(3).rc);
%     c4{subject} = cell2mat(out{subject}{1}.tiss(4).rc);
%     c5{subject} = cell2mat(out{subject}{1}.tiss(5).rc);
%     c6{subject} = cell2mat(out{subject}{1}.tiss(6).rc);
% end

for subject=3:(N+2)
    cd(datadir)
    cd(folders(subject).name);
    
    try
        for tissue_class = 1:6
            tmp = dir(['rc' num2str(tissue_class) '*.nii']);
            if tissue_class == 1
                c1{subject-2} = [tmp.folder filesep tmp.name];
            elseif tissue_class == 2
                c2{subject-2} = [tmp.folder filesep tmp.name];
            elseif tissue_class == 3
                c3{subject-2} = [tmp.folder filesep tmp.name];
            elseif tissue_class == 4
                c4{subject-2} = [tmp.folder filesep tmp.name];
            elseif tissue_class == 5
                c5{subject-2} = [tmp.folder filesep tmp.name];
            elseif tissue_class == 6
                c6{subject-2} = [tmp.folder filesep tmp.name];
            end
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

