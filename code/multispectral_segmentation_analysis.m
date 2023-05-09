%% Data analysis workflow: multispectral segmentation analysis
%
% This is the master script that described and execute the data analysis
% workflow for the 'Model Parametrization for Multispectral Segmentation'
%
% The main aim of the study is to test the SPM12 (r7771) multispectral
% segmentation. To do this we segment data and normalizing
% - using T1 only with 1 Gaussian per tissue class,
% - using T1 only with 2 Gaussians per tissue,
% - using T1 and T2 with 1 Gaussian per tissue class,
% - using T1 and T2 with 2 Gaussians per tissue class
% A Dartel template is also constructed in each case for vizualization.
%
% Several analyses are performed from these images
% - Compare derived volumes of tissue (TIV, GM, WM, CSF and their relarive proportions)
% - Compare tissue distributions via shift funtions and quantify proportions per probability
% - Compare GM images (density using VBM) 

% internal checks 
debug = false; % what is happening: set to another value than false

clear variales
%% set up the directories
root    = fileparts(which('multispectral_segmentation_analysis.m')); cd(root)
datadir = fullfile(root, '..', filesep,'sourcedata', filesep);
outdir  = fullfile(root, '..',  filesep, 'derivatives', filesep);
addpath(root);

export_folder = '';
% '/indirect/staff/cyrilpernet/multispectral_segmentation/Code/SPM_multimodal_segmentation/results/ds003653'
% '/indirect/staff/cyrilpernet/multispectral_segmentation/Code/SPM_multimodal_segmentation/results/NRU_dataset'

if isempty(export_folder)
    error('specify where to export data for statistical analysis and sharing')
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
for op = 1:4
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
    out{op} = segment_images(datadir,outdir,options,debug);
    cd(root); save segmentation_jobs_out out
    
    % compute means and variances
   for class = 1:3
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
    HD{op} = create_decile_images(out{op},outdir,options,debug);
    cd(outdir); save HD HD
    
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



