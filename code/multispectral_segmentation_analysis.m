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
% Add Robust_Statistical_Toolbox-dev to path
addpath(genpath('Robust_Statistical_Toolbox-dev'));

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

%% Compare means and standard deviations between the 4 segmentations
% volumes
load('volumesT1_nG1.mat');  T1_nG1_vol  = volumes; clear volumes
load('volumesT1_nG2.mat');  T1_nG2_vol  = volumes; clear volumes
load('volumesT12_nG1.mat'); T12_nG1_vol = volumes; clear volumes
load('volumesT12_nG2.mat'); T12_nG2_vol = volumes; clear volumes

volumes_std = struct('T1_nG1',std(T1_nG1_vol), 'T1_nG2', std(T1_nG2_vol), ...
                    'T12_nG1', std(T12_nG1_vol), 'T12_nG2', std(T12_nG2_vol));
temp_name = ['volumes' '_std' ];
save(fullfile(outdir, temp_name), 'volumes_std', '-v7.3')

volumes_GM_mean =   [T1_nG1_vol(:,1) T1_nG2_vol(:,1) T12_nG1_vol(:,1) T12_nG2_vol(:,1)];
volumes_mean = struct('T1_nG1',mean(T1_nG1_vol), 'T1_nG2', mean(T1_nG2_vol), ...
                    'T12_nG1', mean(T12_nG1_vol), 'T12_nG2', mean(T12_nG2_vol));
temp_name = ['volumes' '_mean' ];
save(fullfile(outdir, temp_name), 'volumes_mean', '-v7.3')

bar(volumes_mean.T1_nG1);
hold on;
errorbar(volumes_mean.T1_nG1, volumes_std.T1_nG1, '.', 'LineWidth', 1.5);
hold off;
nexttile
bar(volumes_mean.T1_nG2);
hold on;
errorbar(volumes_mean.T1_nG2, volumes_std.T1_nG2, '.', 'LineWidth', 1.5);
hold off;
nexttile
bar(volumes_mean.T12_nG1);
hold on;
errorbar(volumes_mean.T12_nG1, volumes_std.T12_nG1, '.', 'LineWidth', 1.5);
hold off;
nexttile
bar(volumes_mean.T12_nG2);
hold on;
errorbar(volumes_mean.T12_nG2, volumes_std.T12_nG2, '.', 'LineWidth', 1.5);
hold off;

% entropy
load('entropyT1_nG1.mat');  T1_nG1_entropy  = entropy; clear entropy
load('entropyT1_nG2.mat');  T1_nG2_entropy  = entropy; clear entropy
load('entropyT12_nG1.mat'); T12_nG1_entropy = entropy; clear entropy
load('entropyT12_nG2.mat'); T12_nG2_entropy = entropy; clear entropy

entropy_std = struct('T1_nG1',std(T1_nG1_entropy), 'T1_nG2', std(T1_nG2_entropy), ...
                    'T12_nG1', std(T12_nG1_entropy), 'T12_nG2', std(T12_nG2_entropy));
temp_name = ['entropy' '_std' ];
save(fullfile(outdir, temp_name), 'entropy_std', '-v7.3')

entropy_mean = struct('T1_nG1',mean(T1_nG1_entropy), 'T1_nG2', mean(T1_nG2_entropy), ...
                    'T12_nG1', mean(T12_nG1_entropy), 'T12_nG2', mean(T12_nG2_entropy));
temp_name = ['entropy' '_mean' ];
save(fullfile(outdir, temp_name), 'entropy_mean', '-v7.3')

% dunnIndex
load('dunnIndexT1_nG1.mat');  T1_nG1_dunnIndex  = dunnIndexes; clear dunnIndexes
load('dunnIndexT1_nG2.mat');  T1_nG2_dunnIndex  = dunnIndexes; clear dunnIndexes
load('dunnIndexT12_nG1.mat'); T12_nG1_dunnIndex = dunnIndexes; clear dunnIndexes
load('dunnIndexT12_nG2.mat'); T12_nG2_dunnIndex = dunnIndexes; clear dunnIndexes

dunnIndex_std = struct('T1_nG1',std(T1_nG1_dunnIndex), 'T1_nG2', std(T1_nG2_dunnIndex), ...
                    'T12_nG1', std(T12_nG1_dunnIndex), 'T12_nG2', std(T12_nG2_dunnIndex));
temp_name = ['dunnIndex' '_std' ];
save(fullfile(outdir, temp_name), 'dunnIndex_std', '-v7.3')

dunnIndex_mean = struct('T1_nG1',mean(T1_nG1_dunnIndex), 'T1_nG2', mean(T1_nG2_dunnIndex), ...
                    'T12_nG1', mean(T12_nG1_dunnIndex), 'T12_nG2', mean(T12_nG2_dunnIndex));
temp_name = ['dunnIndex' '_mean' ];
save(fullfile(outdir, temp_name), 'dunnIndex_mean', '-v7.3')

% distrib
load('distribT1_nG1.mat');  T1_nG1_distrib  = distrib; clear distrib
load('distribT1_nG2.mat');  T1_nG2_distrib  = distrib; clear distrib
load('distribT12_nG1.mat'); T12_nG1_distrib = distrib; clear distrib
load('distribT12_nG2.mat'); T12_nG2_distrib = distrib; clear distrib

distrib_std = struct('T1_nG1',std(T1_nG1_distrib), 'T1_nG2', std(T1_nG2_distrib), ...
                    'T12_nG1', std(T12_nG1_distrib), 'T12_nG2', std(T12_nG2_distrib));
temp_name = ['distrib' '_std' ];
save(fullfile(outdir, temp_name), 'distrib_std', '-v7.3')

distrib_mean = struct('T1_nG1',mean(T1_nG1_distrib), 'T1_nG2', mean(T1_nG2_distrib), ...
                    'T12_nG1', mean(T12_nG1_distrib), 'T12_nG2', mean(T12_nG2_distrib));
temp_name = ['distrib' '_mean' ];
save(fullfile(outdir, temp_name), 'distrib_mean', '-v7.3')

% distrib_vessels
load('distrib_vesselsT1_nG1.mat');  T1_nG1_distrib_vessels  = distrib_vessels; clear distrib_vessels
load('distrib_vesselsT1_nG2.mat');  T1_nG2_distrib_vessels  = distrib_vessels; clear distrib_vessels
load('distrib_vesselsT12_nG1.mat'); T12_nG1_distrib_vessels = distrib_vessels; clear distrib_vessels
load('distrib_vesselsT12_nG2.mat'); T12_nG2_distrib_vessels = distrib_vessels; clear distrib_vessels

distrib_vessels_std = struct('T1_nG1',std(T1_nG1_distrib_vessels), 'T1_nG2', std(T1_nG2_distrib_vessels), ...
                    'T12_nG1', std(T12_nG1_distrib_vessels), 'T12_nG2', std(T12_nG2_distrib_vessels));
temp_name = ['distrib_vessels' '_std' ];
save(fullfile(outdir, temp_name), 'distrib_vessels_std', '-v7.3')

distrib_vessels_mean = struct('T1_nG1',mean(T1_nG1_distrib_vessels), 'T1_nG2', mean(T1_nG2_distrib_vessels), ...
                    'T12_nG1', mean(T12_nG1_distrib_vessels), 'T12_nG2', mean(T12_nG2_distrib_vessels));
temp_name = ['distrib_vessels' '_mean' ];
save(fullfile(outdir, temp_name), 'distrib_vessels_mean', '-v7.3')

% distrib_nuclei
load('distrib_nucleiT1_nG1.mat');  T1_nG1_distrib_nuclei  = distrib_nuclei; clear distrib_nuclei
load('distrib_nucleiT1_nG2.mat');  T1_nG2_distrib_nuclei  = distrib_nuclei; clear distrib_nuclei
load('distrib_nucleiT12_nG1.mat'); T12_nG1_distrib_nuclei = distrib_nuclei; clear distrib_nuclei
load('distrib_nucleiT12_nG2.mat'); T12_nG2_distrib_nuclei = distrib_nuclei; clear distrib_nuclei

distrib_nuclei_std = struct('T1_nG1',std(T1_nG1_distrib_nuclei), 'T1_nG2', std(T1_nG2_distrib_nuclei), ...
                    'T12_nG1', std(T12_nG1_distrib_nuclei), 'T12_nG2', std(T12_nG2_distrib_nuclei));
temp_name = ['distrib_nuclei' '_std' ];
save(fullfile(outdir, temp_name), 'distrib_nuclei_std', '-v7.3')

distrib_nuclei_mean = struct('T1_nG1',mean(T1_nG1_distrib_nuclei), 'T1_nG2', mean(T1_nG2_distrib_nuclei), ...
                    'T12_nG1', mean(T12_nG1_distrib_nuclei), 'T12_nG2', mean(T12_nG2_distrib_nuclei));
temp_name = ['distrib_nuclei' '_mean' ];
save(fullfile(outdir, temp_name), 'distrib_nuclei_mean', '-v7.3')


%% Compare derived volumes between the 4 segmentations
load('volumesT1_nG1.mat');  T1_nG1  = volumes; clear volumes
load('volumesT1_nG2.mat');  T1_nG2  = volumes; clear volumes
load('volumesT12_nG1.mat'); T12_nG1 = volumes; clear volumes
load('volumesT12_nG2.mat'); T12_nG2 = volumes; clear volumes

%% Compare tissue image distributions between the 4 segmentations
HD
distrib

%% Compare tissue images vozxel-wise between the unimodal and multimodal segmentations with 2 Gaussians
%t-test
