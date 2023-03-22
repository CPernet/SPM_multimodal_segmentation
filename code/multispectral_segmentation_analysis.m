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
