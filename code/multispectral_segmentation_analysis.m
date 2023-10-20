%% Data analysis workflow: multispectral segmentation analysis
%
% This is the master script that describes and executes the data analysis
%
% The main aim of the study is to test the SPM12 (r7771) multispectral
% segmentation. To do this we segment data and normalizing
% - using T1 only with 1 Gaussian per tissue class,
% - using T1 only with 2 Gaussians per tissue,
% - using T1 and T2 with 1 Gaussian per tissue class,
% - using T1 and T2 with 2 Gaussians per tissue class
% A Dartel template is also constructed in each case for vizualization.
%
% A secondary analysis was performed using the new multimodal segmentation
% Multi-Brain https://github.com/WTCN-computational-anatomy-group/mb
%
% Several analyses are performed from these images
% - Compare derived volumes of tissue (TIV, GM, WM, CSF and their relarive proportions)
% - Compare tissue distributions via shift funtions and quantify proportions per probability
% - Compare GM images (density using VBM) 

clear variales

%% set up the directories
root = fileparts(mfilename('fullpath'));
if isempty(root)
    root     = fileparts(which('multispectral_segmentation_analysis.m'));
end
addpath(root); % so we can run the code

% datadir       = fullfile(root(1:strfind(root,'Code')-1),'sourcedata');
% e.g.
datadir       = fullfile(root(1:strfind(root,'Code')-1),'ds003653');
% datadir       = fullfile(root(1:strfind(root,'Code')-1),'nrudataset');

% local output of segmented images and mat files
outdir        = fullfile(datadir, 'derivatives');
if ~exist(outdir,'dir')
    mkdir(outdir)
end
% where we push csv files for sharing
[~,name]=fileparts(datadir);
export_folder = fullfile(fileparts(root),['results' filesep name]);
if ~exist(export_folder,'dir')
    mkdir(export_folder)
end

%% Image processing
% Do the segmentation and get the tissue volumes and voxel distributions. 
% For each case, also make a template image with Dartel (segment_images.m).
% Then create the decile images of tissue density. These relate to the voxel 
% distributions but in space. This is useful to vizualize where are located
% effects observed on distribution of deciles. Finally, mean and variance tissue 
% images are computed allowing a compute a series of Welch's T-tests.

% define the option structure:
% - unimodal vs multimodal: T1 vs T1&T2,
% - number of Gaussians: 1 vs 2
% Jacobian modulation: yes vs no
options = struct('modality',[],'NGaussian',[],'Modulate','No');

% run the segmentation 4 times
for op = 4:-1:1
    if op == 1
        options.modality  = 'T1';
        options.NGaussian = 1;
    elseif op == 2
        options.modality  = 'T1';
        options.NGaussian = 2;
    elseif op == 3
        options.modality  = 'T12';
        options.NGaussian = 1;
    else
        options.modality  = 'T12';
        options.NGaussian = 2;
    end
    out{op} = segment_images(datadir,outdir,options);
    cd(root); % save segmentation_jobs_out out
    
    % compute means and variances
   for class = 3:-1:1
        for n=length(out{op})-1:-1:1
            if op < 3
                V(n) = spm_vol(out{op}{n}{1}.tiss(class).wc{1});
            else
                V(n) = spm_vol(out{op}{n}{2}.tiss(class).wc{1});
            end
        end
        data                = spm_read_vols(V);
        M                   = mean(data,4);
        S                   = var(data,0,4);
        W                   = spm_vol(V(1));
        W.fname             = char(fullfile(outdir,['mean_modality' options.modality '_NGaussian' num2str(options.NGaussian) '_class' num2str(class) '.nii']));
        W.descrip           = 'average image';
        W.private.dat.fname = char(W.fname);
        W.private.dat.dim   = W.dim;
        W.private.descrip   = 'average image';
        W.n                 = [1 1];
        W                   = rmfield(W,'pinfo'); % let SPM figure out the scale
        spm_write_vol(W,M);
        W.fname             = char(fullfile(outdir,['var_modality' options.modality '_NGaussian' num2str(options.NGaussian) '_class' num2str(class) '.nii']));
        W.descrip           = 'variance image';
        W.private.dat.fname = W.fname;
        W.private.descrip   = 'variance image';
        spm_write_vol(W,S);
        clear V data M S W
    end
    
    % get deciles and create decile images
    cd(outdir); 
    if ~exist('HD','var') && exist(fullfile(outdir,'HD.mat'),'file')
        load(fullfile(outdir,'HD.mat')); 
    end
    HD{op} = create_decile_images(out{op},outdir,options);
    save HD HD    
end

%% export most data as csv for sharing/reproducible analyses
cd(outdir)

% volumes
load('volumesT1_nG1.mat');  T1_nG1_vol  = volumes; clear volumes
load('volumesT1_nG2.mat');  T1_nG2_vol  = volumes; clear volumes
load('volumesT12_nG1.mat'); T12_nG1_vol = volumes; clear volumes
load('volumesT12_nG2.mat'); T12_nG2_vol = volumes; clear volumes

volumes_GM  = table(T1_nG1_vol(:,1),T1_nG2_vol(:,1), ...
    T12_nG1_vol(:,1),T12_nG2_vol(:,1), 'VariableNames',...
    {'T1_nG1','T1_nG2','T12_nG1','T12_nG2'});
writetable(volumes_GM,fullfile(export_folder,'GrayMatter_volumes.csv'))

volumes_WM  = table(T1_nG1_vol(:,2),T1_nG2_vol(:,2), ...
    T12_nG1_vol(:,2),T12_nG2_vol(:,2), 'VariableNames',...
    {'T1_nG1','T1_nG2','T12_nG1','T12_nG2'});
writetable(volumes_WM,fullfile(export_folder,'WhiteMatter_volumes.csv'))

volumes_CSF = table(T1_nG1_vol(:,3),T1_nG2_vol(:,3), ...
    T12_nG1_vol(:,3),T12_nG2_vol(:,3), 'VariableNames',...
    {'T1_nG1','T1_nG2','T12_nG1','T12_nG2'});
writetable(volumes_CSF,fullfile(export_folder,'CSF_volumes.csv'))

volumes_Soft = table(T1_nG1_vol(:,4),T1_nG2_vol(:,4), ...
    T12_nG1_vol(:,4),T12_nG2_vol(:,4), 'VariableNames',...
    {'T1_nG1','T1_nG2','T12_nG1','T12_nG2'});
writetable(volumes_Soft,fullfile(export_folder,'SoftTissue_volumes.csv'))

volumes_Skull = table(T1_nG1_vol(:,5),T1_nG2_vol(:,5), ...
    T12_nG1_vol(:,5),T12_nG2_vol(:,5), 'VariableNames',...
    {'T1_nG1','T1_nG2','T12_nG1','T12_nG2'});
writetable(volumes_Skull,fullfile(export_folder,'Skull_volumes.csv'))

volumes_others = table(T1_nG1_vol(:,6),T1_nG2_vol(:,6), ...
    T12_nG1_vol(:,6),T12_nG2_vol(:,6), 'VariableNames',...
    {'T1_nG1','T1_nG2','T12_nG1','T12_nG2'});
writetable(volumes_others,fullfile(export_folder,'Others_volumes.csv'))

% distrib_vessels -- mostly classified as below 0.1 or above 0.9
% summarize across voxels as a frequency
load('distrib_vesselsT1_nG1.mat');  T1_nG1_distrib_vessels  = distrib_vessels; clear distrib_vessels
load('distrib_vesselsT1_nG2.mat');  T1_nG2_distrib_vessels  = distrib_vessels; clear distrib_vessels
load('distrib_vesselsT12_nG1.mat'); T12_nG1_distrib_vessels = distrib_vessels; clear distrib_vessels
load('distrib_vesselsT12_nG2.mat'); T12_nG2_distrib_vessels = distrib_vessels; clear distrib_vessels

tmp1 = T1_nG1_distrib_vessels < 0.1; 
tmp2 = T1_nG2_distrib_vessels < 0.1; 
tmp3 = T12_nG1_distrib_vessels < 0.1;
tmp4 = T12_nG2_distrib_vessels < 0.1;

vessels_GM = table(mean(tmp1(:,:,1),1)'*100,...
    mean(tmp2(:,:,1),1)'*100, mean(tmp3(:,:,1),1)'*100,...
    mean(tmp4(:,:,1),1)'*100, 'VariableNames', ...
    {'T1_nG1','T1_nG2','T12_nG1','T12_nG2'});
writetable(vessels_GM,fullfile(export_folder,'Not_GrayMatter_vessels.csv'))

vessels_WM = table(mean(tmp1(:,:,2),1)'*100,...
    mean(tmp2(:,:,2),1)'*100, mean(tmp3(:,:,2),1)'*100,...
    mean(tmp4(:,:,2),1)'*100, 'VariableNames', ...
    {'T1_nG1','T1_nG2','T12_nG1','T12_nG2'});
writetable(vessels_WM,fullfile(export_folder,'Not_WhiteMatter_vessels.csv'))

vessels_CSF = table(mean(tmp1(:,:,3),1)'*100,...
    mean(tmp2(:,:,3),1)'*100, mean(tmp3(:,:,3),1)'*100,...
    mean(tmp4(:,:,3),1)'*100, 'VariableNames', ...
    {'T1_nG1','T1_nG2','T12_nG1','T12_nG2'});
writetable(vessels_CSF,fullfile(export_folder,'Not_CSF_vessels.csv'))

tmp1 = T1_nG1_distrib_vessels > 0.9; 
tmp2 = T1_nG2_distrib_vessels > 0.9; 
tmp3 = T12_nG1_distrib_vessels > 0.9;
tmp4 = T12_nG2_distrib_vessels > 0.9;

vessels_GM = table(mean(tmp1(:,:,1),1)'*100,...
    mean(tmp2(:,:,1),1)'*100, mean(tmp3(:,:,1),1)'*100,...
    mean(tmp4(:,:,1),1)'*100, 'VariableNames', ...
    {'T1_nG1','T1_nG2','T12_nG1','T12_nG2'});
writetable(vessels_GM,fullfile(export_folder,'GrayMatter_vessels.csv'))

vessels_WM = table(mean(tmp1(:,:,2),1)'*100,...
    mean(tmp2(:,:,2),1)'*100, mean(tmp3(:,:,2),1)'*100,...
    mean(tmp4(:,:,2),1)'*100, 'VariableNames', ...
    {'T1_nG1','T1_nG2','T12_nG1','T12_nG2'});
writetable(vessels_WM,fullfile(export_folder,'WhiteMatter_vessels.csv'))

vessels_CSF = table(mean(tmp1(:,:,3),1)'*100,...
    mean(tmp2(:,:,3),1)'*100, mean(tmp3(:,:,3),1)'*100,...
    mean(tmp4(:,:,3),1)'*100, 'VariableNames', ...
    {'T1_nG1','T1_nG2','T12_nG1','T12_nG2'});
writetable(vessels_CSF,fullfile(export_folder,'CSF_vessels.csv'))

% Dunn Index
load('dunnIndexT1_nG1.mat');  T1_nG1_DI  = dunnIndexes; clear dunnIndexes
load('dunnIndexT1_nG2.mat');  T1_nG2_DI  = dunnIndexes; clear dunnIndexes
load('dunnIndexT12_nG1.mat'); T12_nG1_DI = dunnIndexes; clear dunnIndexes
load('dunnIndexT12_nG2.mat'); T12_nG2_DI = dunnIndexes; clear dunnIndexes

dunnIndexes_GM  = table(T1_nG1_DI(:,1),T1_nG2_DI(:,1), ...
    T12_nG1_DI(:,1),T12_nG2_DI(:,1), 'VariableNames',...
    {'T1_nG1','T1_nG2','T12_nG1','T12_nG2'});
writetable(dunnIndexes_GM,fullfile(export_folder,'GrayMatter_DunnIndexes.csv'))

dunnIndexes_WM  = table(T1_nG1_DI(:,2),T1_nG2_DI(:,2), ...
    T12_nG1_DI(:,2),T12_nG2_DI(:,2), 'VariableNames',...
    {'T1_nG1','T1_nG2','T12_nG1','T12_nG2'});
writetable(dunnIndexes_WM,fullfile(export_folder,'WhiteMatter_DunnIndexes.csv'))

dunnIndexes_CSF = table(T1_nG1_DI(:,3),T1_nG2_DI(:,3), ...
    T12_nG1_DI(:,3),T12_nG2_DI(:,3), 'VariableNames',...
    {'T1_nG1','T1_nG2','T12_nG1','T12_nG2'});
writetable(dunnIndexes_CSF,fullfile(export_folder,'CSF_DunnIndexes.csv'))

% Entropy
load('entropyT1_nG1.mat');  T1_nG1_entropy  = entropy; clear entropy
load('entropyT1_nG2.mat');  T1_nG2_entropy  = entropy; clear entropy
load('entropyT12_nG1.mat'); T12_nG1_entropy = entropy; clear entropy
load('entropyT12_nG2.mat'); T12_nG2_entropy = entropy; clear entropy

entropy_GM  = table(T1_nG1_entropy(:,1),T1_nG2_entropy(:,1), ...
    T12_nG1_entropy(:,1),T12_nG2_entropy(:,1), 'VariableNames',...
    {'T1_nG1','T1_nG2','T12_nG1','T12_nG2'});
writetable(entropy_GM,fullfile(export_folder,'GrayMatter_entropy.csv'))

entropy_WM  = table(T1_nG1_entropy(:,2),T1_nG2_entropy(:,2), ...
    T12_nG1_entropy(:,2),T12_nG2_entropy(:,2), 'VariableNames',...
    {'T1_nG1','T1_nG2','T12_nG1','T12_nG2'});
writetable(entropy_WM,fullfile(export_folder,'WhiteMatter_entropy.csv'))

entropy_CSF = table(T1_nG1_entropy(:,3),T1_nG2_entropy(:,3), ...
    T12_nG1_entropy(:,3),T12_nG2_entropy(:,3), 'VariableNames',...
    {'T1_nG1','T1_nG2','T12_nG1','T12_nG2'});
writetable(entropy_CSF,fullfile(export_folder,'CSF_entropy.csv'))

% just copy HD as it is 
copyfile(fullfile(outdir,'HD.mat'),fullfile(export_folder,'Harrell-Davis-Deciles.mat'))

% just copy nuclei as it is
copyfile(fullfile(outdir,'distrib_nucleiT1_nG1.mat'),fullfile(export_folder,'distrib_nucleiT1_nG1.mat.mat'))
copyfile(fullfile(outdir,'distrib_nucleiT1_nG2.mat'),fullfile(export_folder,'distrib_nucleiT1_nG2.mat.mat'))
copyfile(fullfile(outdir,'distrib_nucleiT12_nG1.mat'),fullfile(export_folder,'distrib_nucleiT12_nG1.mat.mat'))
copyfile(fullfile(outdir,'distrib_nucleiT12_nG2.mat'),fullfile(export_folder,'distrib_nucleiT12_nG2.mat.mat'))
