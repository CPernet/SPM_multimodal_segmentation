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
% Add Robust-Correlation-2 to path
addpath(genpath('Robust-Correlations-2'));

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

%% Compare trimmed means and 95% HDI between the 4 conditions
TissueNames = {'Gray Matter', 'White Matter', 'Cerebrospinal fluid'};
Conditions  = {'T1_nG1', 'T1_nG2', 'T12_nG1', 'T12_nG2', };
% volumes
load('volumesT1_nG1.mat');  T1_nG1_vol  = volumes; clear volumes
load('volumesT1_nG2.mat');  T1_nG2_vol  = volumes; clear volumes
load('volumesT12_nG1.mat'); T12_nG1_vol = volumes; clear volumes
load('volumesT12_nG2.mat'); T12_nG2_vol = volumes; clear volumes

volumes_GM  = [T1_nG1_vol(:,1) T1_nG2_vol(:,1) T12_nG1_vol(:,1) T12_nG2_vol(:,1)];
volumes_WM  = [T1_nG1_vol(:,2) T1_nG2_vol(:,2) T12_nG1_vol(:,2) T12_nG2_vol(:,2)];
volumes_CSF = [T1_nG1_vol(:,3) T1_nG2_vol(:,3) T12_nG1_vol(:,3) T12_nG2_vol(:,3)];

[GM_est, CI_GM]   = rst_data_plot(volumes_GM, 'estimator','trimmed mean');
[WM_est, CI_WM]   = rst_data_plot(volumes_WM, 'estimator','trimmed mean','newfig','yes');
[CSF_est, CI_CSF] = rst_data_plot(volumes_CSF, 'estimator','trimmed mean','newfig','yes');

% Add x-axis labels
XL       = get(findall(figure(1),'type','axes'), 'XLim');
ticLengh = ((XL(2)-XL(1))/4);
xticks(ticLengh-XL(1) : ticLengh : (ticLengh*4)-XL(1));
xticklabels(Conditions);
XL       = get(findall(figure(2),'type','axes'), 'XLim');
ticLengh = ((XL(2)-XL(1))/4);
xticks(ticLengh-XL(1) : ticLengh : (ticLengh*4)-XL(1));
xticklabels(Conditions);
XL       = get(findall(figure(3),'type','axes'), 'XLim');
ticLengh = ((XL(2)-XL(1))/4);
xticks(ticLengh-XL(1) : ticLengh : (ticLengh*4)-XL(1));
xticklabels(Conditions);
if(debug)
    % save plot for volumes trimmed mean and close figure
    saveas(figure(1), "GM_Volumes_TM.png"); 
    saveas(figure(2), "WM_Volumes_TM.png");
    saveas(figure(3), "CSF_Volumes_TM.png");
end
close(figure(1));
close(figure(2));
close(figure(3));

TrimmedMeans = [GM_est; WM_est; CSF_est];
LowerConfs   = [CI_GM(1,:); CI_WM(1,:); CI_CSF(1,:)];
HigherConfs  = [CI_GM(2,:); CI_WM(2,:); CI_CSF(2,:)];

T1_nG1  = [LowerConfs(:,1) TrimmedMeans(:,1) HigherConfs(:,1)];
T1_nG2  = [LowerConfs(:,2) TrimmedMeans(:,2) HigherConfs(:,2)];
T12_nG1 = [LowerConfs(:,3) TrimmedMeans(:,3) HigherConfs(:,3)];
T12_nG2 = [LowerConfs(:,4) TrimmedMeans(:,4) HigherConfs(:,4)];

T = table(T1_nG1,T1_nG2,T12_nG1,T12_nG2,...
    'RowNames',TissueNames);
writetable(T,'Volumes_TrimmedMeans.csv','WriteRowNames',true);
clear GM_est WM_est CSF_est CI_GM CI_WM CI_CSF TrimmedMeans LowerConfs HigherConfs T1_nG1 T1_nG2 T12_nG1 T12_nG2 T XL ticLengh

% Multi compare between the 4 conditions (T1_nG1 vs T1_nG2, T12_nG1 vs T12_nG2, T1_nG1 vs T12_nG1, T1_nG2 vs T12_nG2)
[diff_GM,CI_GM,p_GM,alphav_GM,h_GM]      = rst_multicompare(volumes_GM,[1 2; 3 4; 1 3; 2 4], 'estimator', 'trimmed mean','newfig','yes');
[diff_WM,CI_WM,p_WM,alphav_WM,h_WM]      = rst_multicompare(volumes_WM,[1 2; 3 4; 1 3; 2 4], 'estimator', 'trimmed mean','newfig','yes');
[diff_CSF,CI_CSF,p_CSF,alphav_CSF,h_CSF] = rst_multicompare(volumes_CSF,[1 2; 3 4; 1 3; 2 4], 'estimator', 'trimmed mean','newfig','yes');
% Add x-axis labels
XL       = get(findall(figure(1),'type','axes'), 'XLim');
ticLengh = ((XL(2)-XL(1))/4);
xticks(ticLengh : ticLengh : (ticLengh*4));
xticklabels(Conditions);
XL       = get(findall(figure(2),'type','axes'), 'XLim');
ticLengh = ((XL(2)-XL(1))/4);
xticks(ticLengh : ticLengh : (ticLengh*4));
xticklabels(Conditions);
XL       = get(findall(figure(3),'type','axes'), 'XLim');
ticLengh = ((XL(2)-XL(1))/4);
xticks(ticLengh : ticLengh : (ticLengh*4));
xticklabels(Conditions);
if(debug)
    % save plot for Multi compare and close figure
    saveas(figure(1), "MultiComp_GM.png");
    saveas(figure(2), "MultiComp_WM.png");
    saveas(figure(3), "MultiComp_CSF.png");
end
close(figure(1));
close(figure(2));
close(figure(3));

PairwiseDifferences = [diff_GM; diff_WM; diff_CSF];
LowerConfs          = [CI_GM(1,:); CI_WM(1,:); CI_CSF(1,:)];
HigherConfs         = [CI_GM(2,:); CI_WM(2,:); CI_CSF(2,:)];
PValues             = [p_GM; p_WM; p_CSF];
AlphaValues         = [alphav_GM; alphav_WM; alphav_CSF];
Significances       = [h_GM'; h_WM'; h_CSF'];

T = table(PairwiseDifferences,LowerConfs,HigherConfs, ...
          PValues,AlphaValues,Significances,...
          'RowNames',TissueNames);
writetable(T,'Volumes_multicompare.csv','WriteRowNames',true);
clear diff_GM diff_WM diff_CSF CI_GM CI_WM CI_CSF p_GM p_WM p_CSF alphav_GM alphav_WM h_GM h_WM h_CSF alphav_CSF PairwiseDifferences LowerConfs HigherConfs PValues AlphaValues Significances T

% correlations GM vs WM GM vs CSF WM vs CSF
[rp_T1_nG1,tp_T1_nG1,CI_T1_nG1,pval_T1_nG1,outid_T1_nG1,h_T1_nG1]       = skipped_Spearman(T1_nG1_vol,[1 2; 1 3; 3 2]);
[rp_T1_nG2,tp_T1_nG2,CI_T1_nG2,pval_T1_nG2,outid_T1_nG2,h_T1_nG2]       = skipped_Spearman(T1_nG2_vol,[1 2; 1 3; 3 2]);
[rp_T12_nG1,tp_T12_nG1,CI_T12_nG1,pval_T12_nG1,outid_T12_nG1,h_T12_nG1] = skipped_Spearman(T12_nG1_vol,[1 2; 1 3; 3 2]);
[rp_T12_nG2,tp_T12_nG2,CI_T12_nG2,pval_T12_nG2,outid_T12_nG2,h_T12_nG2] = skipped_Spearman(T12_nG2_vol,[1 2; 1 3; 3 2]);

SpearmanCorrelations = [rp_T1_nG1'; rp_T1_nG2'; rp_T12_nG1'; rp_T12_nG2'];
TValues              = [tp_T1_nG1'; tp_T1_nG2'; tp_T12_nG1'; tp_T12_nG2'];
LowerConfs           = [CI_T1_nG1(1,:); CI_T1_nG2(1,:); CI_T12_nG1(1,:); CI_T12_nG2(1,:)];
HigherConfs          = [CI_T1_nG1(2,:); CI_T1_nG2(2,:); CI_T12_nG1(2,:); CI_T12_nG2(2,:)];
PValues              = [pval_T1_nG1'; pval_T1_nG2'; pval_T12_nG1'; pval_T12_nG2'];
IndexOutliers        = [outid_T1_nG1'; outid_T1_nG2'; outid_T12_nG1'; outid_T12_nG2'];
Significances        = [h_T1_nG1'; h_T1_nG2'; h_T12_nG1'; h_T12_nG2'];

T = table(SpearmanCorrelations,TValues,LowerConfs, ...
          HigherConfs,PValues,IndexOutliers, Significances,...
          'RowNames',Conditions);
writetable(T,'Volumes_spearman.csv','WriteRowNames',true);
clear rp_T1_nG1 tp_T1_nG1 CI_T1_nG1 pval_T1_nG1 outid_T1_nG1 h_T1_nG1 rp_T1_nG2 tp_T1_nG2 CI_T1_nG2 pval_T1_nG2 outid_T1_nG2 h_T1_nG2 rp_T12_nG1 tp_T12_nG1 CI_T12_nG1 pval_T12_nG1 outid_T12_nG1 h_T12_nG1 rp_T12_nG2 tp_T12_nG2 CI_T12_nG2 pval_T12_nG2 outid_T12_nG2 h_T12_nG2 SpearmanCorrelations TValues LowerConfs HigherConfs PValues IndexOutliers Significances T TissueNames ticLengh XL

% Total IOntracranial Volume (GM+WM+CSF) -- link to dartel template, where
% does extra-missing volumes go?
TIV = volumes_GM + volumes_WM + volumes_CSF;
[T1_nG1_est, CI_T1_nG1]   = rst_data_plot(TIV(:,1), 'estimator','trimmed mean');
[T1_nG2_est, CI_T1_nG2]   = rst_data_plot(TIV(:,2), 'estimator','trimmed mean','newfig','yes');
[T12_nG1_est, CI_T12_nG1] = rst_data_plot(TIV(:,3), 'estimator','trimmed mean','newfig','yes');
[T12_nG2_est, CI_T12_nG2] = rst_data_plot(TIV(:,4), 'estimator','trimmed mean','newfig','yes');
% Add x-axis labels
XL       = get(findall(figure(1),'type','axes'), 'XLim');
ticLengh = ((XL(2))/2);
xticks(ticLengh);
xticklabels(Conditions{1});
XL       = get(findall(figure(2),'type','axes'), 'XLim');
ticLengh = ((XL(2))/2);
xticks(ticLengh);
xticklabels(Conditions{2});
XL       = get(findall(figure(3),'type','axes'), 'XLim');
ticLengh = ((XL(2))/2);
xticks(ticLengh);
xticklabels(Conditions{3});
XL       = get(findall(figure(4),'type','axes'), 'XLim');
ticLengh = ((XL(2))/2);
xticks(ticLengh);
xticklabels(Conditions{4});
if(debug)
    % save plot for volumes trimmed mean and close figure
    saveas(figure(1), "TIV_T1_nG1_Volume_TM.png"); 
    saveas(figure(2), "TIV_T1_nG2_Volume_TM.png");
    saveas(figure(3), "TIV_T12_nG1_Volume_TM.png"); 
    saveas(figure(4), "TIV_T12_nG2_Volume_TM.png");
end
close(figure(1));
close(figure(2));
close(figure(3));
close(figure(4));

TrimmedMeans = [T1_nG1_est, T1_nG2_est, T12_nG1_est, T12_nG2_est];
LowerConfs   = [CI_T1_nG1(1,:), CI_T1_nG2(1,:), CI_T12_nG1(1,:), CI_T12_nG2(1,:)];
HigherConfs  = [CI_T1_nG1(2,:), CI_T1_nG2(2,:), CI_T12_nG1(2,:), CI_T12_nG2(2,:)];

T1_nG1  = [LowerConfs(:,1) TrimmedMeans(:,1) HigherConfs(:,1)];
T1_nG2  = [LowerConfs(:,2) TrimmedMeans(:,2) HigherConfs(:,2)];
T12_nG1 = [LowerConfs(:,3) TrimmedMeans(:,3) HigherConfs(:,3)];
T12_nG2 = [LowerConfs(:,4) TrimmedMeans(:,4) HigherConfs(:,4)];

T = table(T1_nG1,T1_nG2,T12_nG1,T12_nG2);
writetable(T,'TIV_Volume_TrimmedMeans.csv','WriteRowNames',true);
clear T1_nG1_est CI_T1_nG1 T1_nG2_est CI_T1_nG2 T12_nG1_est CI_T12_nG1 T12_nG2_est CI_T12_nG2 TrimmedMeans LowerConfs HigherConfs T1_nG1 T1_nG2 T12_nG1 T12_nG2 T volumes_GM volumes_WM volumes_CSF T1_nG1_vol T1_nG2_vol T12_nG1_vol T12_nG2_vol XL ticLengh TIV

% distrib_vessels
load('distrib_vesselsT1_nG1.mat');  T1_nG1_distrib_vessels  = distrib_vessels; clear distrib_vessels
load('distrib_vesselsT1_nG2.mat');  T1_nG2_distrib_vessels  = distrib_vessels; clear distrib_vessels
load('distrib_vesselsT12_nG1.mat'); T12_nG1_distrib_vessels = distrib_vessels; clear distrib_vessels
load('distrib_vesselsT12_nG2.mat'); T12_nG2_distrib_vessels = distrib_vessels; clear distrib_vessels

% HD
load('HD.mat');

% entropy
load('entropyT1_nG1.mat');  T1_nG1_entropy  = entropy; clear entropy
load('entropyT1_nG2.mat');  T1_nG2_entropy  = entropy; clear entropy
load('entropyT12_nG1.mat'); T12_nG1_entropy = entropy; clear entropy
load('entropyT12_nG2.mat'); T12_nG2_entropy = entropy; clear entropy

entropy_GM  = [T1_nG1_entropy(:,1) T1_nG2_entropy(:,1) T12_nG1_entropy(:,1) T12_nG2_entropy(:,1)];
entropy_WM  = [T1_nG1_entropy(:,2) T1_nG2_entropy(:,2) T12_nG1_entropy(:,2) T12_nG2_entropy(:,2)];
entropy_CSF = [T1_nG1_entropy(:,3) T1_nG2_entropy(:,3) T12_nG1_entropy(:,3) T12_nG2_entropy(:,3)];

[GM_est, CI_GM]   = rst_data_plot(entropy_GM, 'estimator','trimmed mean');
[WM_est, CI_WM]   = rst_data_plot(entropy_WM, 'estimator','trimmed mean','newfig','yes');
[CSF_est, CI_CSF] = rst_data_plot(entropy_CSF, 'estimator','trimmed mean','newfig','yes');
if(debug)
    % save plot for entropy trimmed mean and close figure
    saveas(figure(1), "GM_entropy_TM.png"); 
    saveas(figure(2), "WM_entropy_TM.png");
    saveas(figure(3), "CSF_entropy_TM.png");
end
close(figure(1));
close(figure(2));
close(figure(3));

TrimmedMeans = [GM_est; WM_est; CSF_est];
LowerConfs   = [CI_GM(1,:); CI_WM(1,:); CI_CSF(1,:)];
HigherConfs  = [CI_GM(2,:); CI_WM(2,:); CI_CSF(2,:)];

T1_nG1  = [LowerConfs(:,1) TrimmedMeans(:,1) HigherConfs(:,1)];
T1_nG2  = [LowerConfs(:,2) TrimmedMeans(:,2) HigherConfs(:,2)];
T12_nG1 = [LowerConfs(:,3) TrimmedMeans(:,3) HigherConfs(:,3)];
T12_nG2 = [LowerConfs(:,4) TrimmedMeans(:,4) HigherConfs(:,4)];

T = table(T1_nG1,T1_nG2,T12_nG1,T12_nG2,...
    'RowNames',RowNames);
writetable(T,'entropy_TrimmedMeans.csv','WriteRowNames',true);
clear GM_est WM_est CSF_est CI_GM CI_WM CI_CSF TrimmedMeans LowerConfs HigherConfs T1_nG1 T1_nG2 T12_nG1 T12_nG2 T

% Multi compare between the 4 conditions (T1_nG1 vs T1_nG2, T12_nG1 vs T12_nG2, T1_nG1 vs T12_nG1, T1_nG2 vs T12_nG2)
[diff_GM,CI_GM,p_GM,alphav_GM,h_GM]      = rst_multicompare(entropy_GM,[1 2; 3 4; 1 3; 2 4], 'estimator', 'trimmed mean','newfig','yes');
[diff_WM,CI_WM,p_WM,alphav_WM,h_WM]      = rst_multicompare(entropy_WM,[1 2; 3 4; 1 3; 2 4], 'estimator', 'trimmed mean','newfig','yes');
[diff_CSF,CI_CSF,p_CSF,alphav_CSF,h_CSF] = rst_multicompare(entropy_CSF,[1 2; 3 4; 1 3; 2 4], 'estimator', 'trimmed mean','newfig','yes');
if(debug)
    % save plot for Multi compare and close figure
    saveas(figure(1), "MultiComp_GM_entropy.png");
    saveas(figure(2), "MultiComp_WM_entropy.png");
    saveas(figure(3), "MultiComp_CSF_entropy.png");
end
close(figure(1));
close(figure(2));
close(figure(3));

PairwiseDifferences = [diff_GM; diff_WM; diff_CSF];
LowerConfs          = [CI_GM(1,:); CI_WM(1,:); CI_CSF(1,:)];
HigherConfs         = [CI_GM(2,:); CI_WM(2,:); CI_CSF(2,:)];
PValues             = [p_GM,p_WM, p_CSF];
AlphaValues         = [alphav_GM; alphav_WM; alphav_CSF];
Significances       = [h_GM; h_WM; h_CSF];

clear diff_GM diff_WM diff_CSF CI_GM CI_WM CI_CSF PairwiseDifferences LowerConfs HigherConfs PValues AlphaValues Significances T


% dunnIndex
load('dunnIndexT1_nG1.mat');  T1_nG1_dunnIndex  = dunnIndexes; clear dunnIndexes
load('dunnIndexT1_nG2.mat');  T1_nG2_dunnIndex  = dunnIndexes; clear dunnIndexes
load('dunnIndexT12_nG1.mat'); T12_nG1_dunnIndex = dunnIndexes; clear dunnIndexes
load('dunnIndexT12_nG2.mat'); T12_nG2_dunnIndex = dunnIndexes; clear dunnIndexes

dunnIndex_GM  = [T1_nG1_dunnIndex(:,1) T1_nG2_dunnIndex(:,1) T12_nG1_dunnIndex(:,1) T12_nG2_dunnIndex(:,1)];
dunnIndex_WM  = [T1_nG1_dunnIndex(:,2) T1_nG2_dunnIndex(:,2) T12_nG1_dunnIndex(:,2) T12_nG2_dunnIndex(:,2)];
dunnIndex_CSF = [T1_nG1_dunnIndex(:,3) T1_nG2_dunnIndex(:,3) T12_nG1_dunnIndex(:,3) T12_nG2_dunnIndex(:,3)];

[GM_est, CI_GM]   = rst_data_plot(dunnIndex_GM, 'estimator','trimmed mean');
[WM_est, CI_WM]   = rst_data_plot(dunnIndex_WM, 'estimator','trimmed mean','newfig','yes');
[CSF_est, CI_CSF] = rst_data_plot(dunnIndex_CSF, 'estimator','trimmed mean','newfig','yes');
if(debug)
    % save plot for dunnIndex trimmed mean and close figure
    saveas(figure(1), "GM_dunnIndex_TM.png"); 
    saveas(figure(2), "WM_dunnIndex_TM.png");
    saveas(figure(3), "CSF_dunnIndex_TM.png");
end
close(figure(1));
close(figure(2));
close(figure(3));

TrimmedMeans = [GM_est; WM_est; CSF_est];
LowerConfs   = [CI_GM(1,:); CI_WM(1,:); CI_CSF(1,:)];
HigherConfs  = [CI_GM(2,:); CI_WM(2,:); CI_CSF(2,:)];

T1_nG1  = [LowerConfs(:,1) TrimmedMeans(:,1) HigherConfs(:,1)];
T1_nG2  = [LowerConfs(:,2) TrimmedMeans(:,2) HigherConfs(:,2)];
T12_nG1 = [LowerConfs(:,3) TrimmedMeans(:,3) HigherConfs(:,3)];
T12_nG2 = [LowerConfs(:,4) TrimmedMeans(:,4) HigherConfs(:,4)];

T = table(T1_nG1,T1_nG2,T12_nG1,T12_nG2,...
    'RowNames',RowNames);
writetable(T,'dunnIndex_TrimmedMeans.csv','WriteRowNames',true);
clear GM_est WM_est CSF_est CI_GM CI_WM CI_CSF TrimmedMeans LowerConfs HigherConfs T1_nG1 T1_nG2 T12_nG1 T12_nG2 T

% Multi compare between the 4 conditions (T1_nG1 vs T1_nG2, T12_nG1 vs T12_nG2, T1_nG1 vs T12_nG1, T1_nG2 vs T12_nG2)
[diff_GM,CI_GM,p_GM,alphav_GM,h_GM]      = rst_multicompare(dunnIndex_GM,[1 2; 3 4; 1 3; 2 4], 'estimator', 'trimmed mean','newfig','yes');
[diff_WM,CI_WM,p_WM,alphav_WM,h_WM]      = rst_multicompare(dunnIndex_WM,[1 2; 3 4; 1 3; 2 4], 'estimator', 'trimmed mean','newfig','yes');
[diff_CSF,CI_CSF,p_CSF,alphav_CSF,h_CSF] = rst_multicompare(dunnIndex_CSF,[1 2; 3 4; 1 3; 2 4], 'estimator', 'trimmed mean','newfig','yes');
if(debug)
    % save plot for Multi compare and close figure
    saveas(figure(1), "MultiComp_GM_dunnIndex.png");
    saveas(figure(2), "MultiComp_WM_dunnIndex.png");
    saveas(figure(3), "MultiComp_CSF_dunnIndex.png");
end
close(figure(1));
close(figure(2));
close(figure(3));

PairwiseDifferences = [diff_GM; diff_WM; diff_CSF];
LowerConfs          = [CI_GM(1,:); CI_WM(1,:); CI_CSF(1,:)];
HigherConfs         = [CI_GM(2,:); CI_WM(2,:); CI_CSF(2,:)];
PValues             = [p_GM,p_WM, p_CSF];
AlphaValues         = [alphav_GM; alphav_WM; alphav_CSF];
Significances       = [h_GM; h_WM; h_CSF];

clear diff_GM diff_WM diff_CSF CI_GM CI_WM CI_CSF PairwiseDifferences LowerConfs HigherConfs PValues AlphaValues Significances T




% distrib
load('distribT1_nG1.mat');  T1_nG1_distrib  = distrib; clear distrib
load('distribT1_nG2.mat');  T1_nG2_distrib  = distrib; clear distrib
load('distribT12_nG1.mat'); T12_nG1_distrib = distrib; clear distrib
load('distribT12_nG2.mat'); T12_nG2_distrib = distrib; clear distrib

tissue_distrib_GM  = [T1_nG1_distrib(:,1) T1_nG2_distrib(:,1) T12_nG1_distrib(:,1) T12_nG2_distrib(:,1)];
tissue_distrib_WM  = [T1_nG1_distrib(:,2) T1_nG2_distrib(:,2) T12_nG1_distrib(:,2) T12_nG2_distrib(:,2)];
tissue_distrib_CSF = [T1_nG1_distrib(:,3) T1_nG2_distrib(:,3) T12_nG1_distrib(:,3) T12_nG2_distrib(:,3)];

[GM_est, CI_GM]   = rst_data_plot(tissue_distrib_GM, 'estimator','trimmed mean');
[WM_est, CI_WM]   = rst_data_plot(tissue_distrib_WM, 'estimator','trimmed mean','newfig','yes');
[CSF_est, CI_CSF] = rst_data_plot(tissue_distrib_CSF, 'estimator','trimmed mean','newfig','yes');
if(debug)
    % save plot for tissue distrib trimmed mean and close figure
    saveas(figure(1), "GM_tissue_distrib_TM.png"); 
    saveas(figure(2), "WM_tissue_distrib_TM.png");
    saveas(figure(3), "CSF_tissue_distrib_TM.png");
end
close(figure(1));
close(figure(2));
close(figure(3));

TrimmedMeans = [GM_est; WM_est; CSF_est];
LowerConfs   = [CI_GM(1,:); CI_WM(1,:); CI_CSF(1,:)];
HigherConfs  = [CI_GM(2,:); CI_WM(2,:); CI_CSF(2,:)];

T1_nG1  = [LowerConfs(:,1) TrimmedMeans(:,1) HigherConfs(:,1)];
T1_nG2  = [LowerConfs(:,2) TrimmedMeans(:,2) HigherConfs(:,2)];
T12_nG1 = [LowerConfs(:,3) TrimmedMeans(:,3) HigherConfs(:,3)];
T12_nG2 = [LowerConfs(:,4) TrimmedMeans(:,4) HigherConfs(:,4)];

T = table(T1_nG1,T1_nG2,T12_nG1,T12_nG2,...
    'RowNames',RowNames);
writetable(T,'tissue_distrib_TrimmedMeans.csv','WriteRowNames',true);
clear GM_est WM_est CSF_est CI_GM CI_WM CI_CSF TrimmedMeans LowerConfs HigherConfs T1_nG1 T1_nG2 T12_nG1 T12_nG2 T

% Multi compare between the 4 conditions (T1_nG1 vs T1_nG2, T12_nG1 vs T12_nG2, T1_nG1 vs T12_nG1, T1_nG2 vs T12_nG2)
[diff_GM,CI_GM,p_GM,alphav_GM,h_GM] = rst_multicompare(tissue_distrib_GM,[1 2; 3 4; 1 3; 2 4], 'estimator', 'trimmed mean','newfig','yes');
[diff_WM,CI_WM,p_WM,alphav_WM,h_WM] = rst_multicompare(tissue_distrib_WM,[1 2; 3 4; 1 3; 2 4], 'estimator', 'trimmed mean','newfig','yes');
[diff_CSF,CI_CSF,p_CSF,alphav_CSF,h_CSF] = rst_multicompare(tissue_distrib_CSF,[1 2; 3 4; 1 3; 2 4], 'estimator', 'trimmed mean','newfig','yes');
if(debug)
    % save plot for Multi compare and close figure
    saveas(figure(1), "MultiComp_GM_tissue_distrib.png");
    saveas(figure(2), "MultiComp_WM_tissue_distrib.png");
    saveas(figure(3), "MultiComp_CSF_tissue_distrib.png");
end
close(figure(1));
close(figure(2));
close(figure(3));

PairwiseDifferences = [diff_GM; diff_WM; diff_CSF];
LowerConfs          = [CI_GM(1,:); CI_WM(1,:); CI_CSF(1,:)];
HigherConfs         = [CI_GM(2,:); CI_WM(2,:); CI_CSF(2,:)];
PValues             = [p_GM,p_WM, p_CSF];
AlphaValues         = [alphav_GM; alphav_WM; alphav_CSF];
Significances       = [h_GM; h_WM; h_CSF];

clear diff_GM diff_WM diff_CSF CI_GM CI_WM CI_CSF PairwiseDifferences LowerConfs HigherConfs PValues AlphaValues Significances T




% distrib_nuclei
load('distrib_nucleiT1_nG1.mat');  T1_nG1_distrib_nuclei  = distrib_nuclei; clear distrib_nuclei
load('distrib_nucleiT1_nG2.mat');  T1_nG2_distrib_nuclei  = distrib_nuclei; clear distrib_nuclei
load('distrib_nucleiT12_nG1.mat'); T12_nG1_distrib_nuclei = distrib_nuclei; clear distrib_nuclei
load('distrib_nucleiT12_nG2.mat'); T12_nG2_distrib_nuclei = distrib_nuclei; clear distrib_nuclei


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
