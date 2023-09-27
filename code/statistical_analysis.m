%% Statistical Analysis of multispectral segmentation 
% 
% The analysis proceed such as NRU data (N=259) are the discovery set, 
% while the public OpenNeuro ds003653 data (N=87) are the validation set.
% For each analysis, we compute effect sizes with confidence intervals to
% evaluate if effects can reproduce and test for significance using NHST to
% evaluate if effects can replicate.

addpath('external')
cd('../results')

% read the data
GMd  = readtable(['nrudataset' filesep 'GrayMatter_volumes.csv'],'ReadRowNames',false);           
WMd  = readtable(['nrudataset' filesep 'WhiteMatter_volumes.csv'],'ReadRowNames',false);           
CSFd = readtable(['nrudataset' filesep 'CSF_volumes.csv'],'ReadRowNames',false);           
GMt  = readtable(['ds003653' filesep 'GrayMatter_volumes.csv'],'ReadRowNames',false);           
WMt  = readtable(['ds003653' filesep 'WhiteMatter_volumes.csv'],'ReadRowNames',false);           
CSFt = readtable(['ds003653' filesep 'CSF_volumes.csv'],'ReadRowNames',false);           

%% Volume analyses

% ----------------------------------------- 
%% 1 - volume estimates and reproducibility
% ------------------------------------------

TIVd = [GMd{:,1}+WMd{:,1}+CSFd{:,1} GMd{:,2}+WMd{:,2}+CSFd{:,2} ...
    GMd{:,3}+WMd{:,3}+CSFd{:,3} GMd{:,4}+WMd{:,4}+CSFd{:,4}].*1000;
TIVt = [GMt{:,1}+WMt{:,1}+CSFt{:,1} GMt{:,2}+WMt{:,2}+CSFt{:,2} ...
    GMt{:,3}+WMt{:,3}+CSFt{:,3} GMt{:,4}+WMt{:,4}+CSFt{:,4}].*1000;

% start figure for TIV and get estimates
figure('Name','Tissue volumes'); subplot(6,13,[1 2 3 13+1 13+2 13+3 26+1 26+2 26+3]);
[TIVd_est, TIVd_CI]   = rst_data_plot(TIVd, 'estimator','trimmed mean','newfig','sub');
title('TIV discovery set','Fontsize',12); ylabel('volumes'); axis([0.1 5.9 990 1890])
subplot(6,13,[39+1 39+2 39+3 52+1 52+2 52+3 65+1 65+2 65+3]);
[TIVt_est, TIVt_CI]   = rst_data_plot(TIVt, 'estimator','trimmed mean','newfig','sub');
title('TIV validation set','Fontsize',12); ylabel('volumes'); axis([0.1 5.9 990 1890])
xlabel('T1-1G T1-2G T12-1G T12-2G'); 

% get estimates for GM, WM and CSF, make a supplementary figure
figure('Name','Tissue Class volumes'); subplot(2,3,1);
[GMd_est, CId_GM,~,K1]   = rst_data_plot(GMd{:,:}.*1000, 'estimator','trimmed mean','newfig','sub');
title('Grey Matter discovery set','Fontsize',12); ylabel('GM volumes'); axis([0.1 5.9 490 990]); subplot(2,3,2);
[WMd_est, CId_WM,~,K2]   = rst_data_plot(WMd{:,:}.*1000, 'estimator','trimmed mean','newfig','sub');
title('White Matter discovery set','Fontsize',12); ylabel('WM volumes'); axis([0.1 5.9 290 690]); subplot(2,3,3);
[CSFd_est, CId_CSF,~,K3] = rst_data_plot(CSFd{:,:}.*1000, 'estimator','trimmed mean','newfig','sub');
title('CSF discovery set','Fontsize',12); ylabel('CSF volumes'); axis([0.1 5.9 90 430]); subplot(2,3,4);
[GMt_est, CIt_GM,~,K4]   = rst_data_plot(GMt{:,:}.*1000, 'estimator','trimmed mean','newfig','sub');
title('Grey Matter test set','Fontsize',12); ylabel('GM volumes'); axis([0.1 5.9 490 990]); subplot(2,3,5);
[WMt_est, CIt_WM,~,K5]   = rst_data_plot(WMt{:,:}.*1000, 'estimator','trimmed mean','newfig','sub');
title('White Matter test set','Fontsize',12); ylabel('WM volumes'); axis([0.1 5.9 290 690]); subplot(2,3,6);
[CSFt_est, CIt_CSF,~,K6] = rst_data_plot(CSFt{:,:}.*1000, 'estimator','trimmed mean','newfig','sub');
title('CSF test set','Fontsize',12); ylabel('CSF volumes'); axis([0.1 5.9 90 430]);

% complete figure 1 using scatter plots and kernel density estimates
figure(findobj( 'Type', 'Figure', 'Name', 'Tissue volumes' ));

% discovery
gp = [repmat({'1 Gaussian'},size(GMd,1),1);repmat({'2 Gaussians'},size(GMd,1),1)]; 
subplot(6,13,[5 6]); plot(K1{1},'b','LineWidth',2); hold on; plot(K1{2},'r','LineWidth',2); 
grid on; set(gca,'xticklabel',{[]}); title('GM volumes')
subplot(6,13,[13+5 13+6 26+5 26+6]); gscatter([GMd{:,1};GMd{:,2}].*1000,[GMd{:,3};GMd{:,4}].*1000,gp);
grid on; axis('square'); xlabel('T1w'); ylabel('T1w & T2w'); legend('off'); axis([400 1000 400 1000])
subplot(6,13,[13+7 26+7]); plot(fliplr(K1{3}),'b','LineWidth',2); hold on; 
grid on; plot(fliplr(K1{4}),'r','LineWidth',2); set(gca,'xticklabel',{[]}); camroll(-90)

subplot(6,13,[8 9]); plot(K2{1},'b','LineWidth',2); hold on; plot(K2{2},'r','LineWidth',2); 
grid on; set(gca,'xticklabel',{[]}); title('WM volumes')
subplot(6,13,[13+8 13+9 26+8 26+9]); gscatter([WMd{:,1};WMd{:,2}].*1000,[WMd{:,3};WMd{:,4}].*1000,gp);
grid on; axis('square'); xlabel('T1w'); legend('off'); axis([300 650 300 650])
subplot(6,13,[13+10 26+10]); plot(fliplr(K2{3}),'b','LineWidth',2); hold on; 
grid on; plot(fliplr(K2{4}),'r','LineWidth',2); set(gca,'xticklabel',{[]}); camroll(-90)

subplot(6,13,[11 12]); plot(K3{1},'b','LineWidth',2); hold on; plot(K3{2},'r','LineWidth',2); 
grid on; set(gca,'xticklabel',{[]}); title('CSF volumes')
subplot(6,13,[13+11 13+12 26+11 26+12]); gscatter([CSFd{:,1};CSFd{:,2}].*1000,[CSFd{:,3};CSFd{:,4}].*1000,gp);
grid on; axis('square'); xlabel('T1w'); legend('off'); axis([100 450 100 450])
subplot(6,13,[13+13 26+13]); plot(fliplr(K3{3}),'b','LineWidth',2); hold on; 
grid on; plot(fliplr(K3{4}),'r','LineWidth',2); set(gca,'xticklabel',{[]}); camroll(-90)

% validatation
gp = [repmat({'1 Gaussian'},size(GMt,1),1);repmat({'2 Gaussians'},size(GMt,1),1)]; 
subplot(6,13,[39+5 39+6]); plot(K4{1},'b','LineWidth',2); hold on; plot(K4{2},'r','LineWidth',2); 
grid on; set(gca,'xticklabel',{[]}); 
subplot(6,13,[52+5 52+6 65+5 65+6]); gscatter([GMt{:,1};GMt{:,2}].*1000,[GMt{:,3};GMt{:,4}].*1000,gp);
grid on; axis('square'); xlabel('T1w'); ylabel('T1w & T2w'); legend('off'); axis([400 1000 400 1000])
subplot(6,13,[52+7 65+7]); plot(K1{3},'b','LineWidth',2); hold on; 
grid on; plot(K1{4},'r','LineWidth',2); set(gca,'xticklabel',{[]}); camroll(-90)

subplot(6,13,[39+8 39+9]); plot(K5{1},'b','LineWidth',2); hold on; plot(K5{2},'r','LineWidth',2); 
grid on; set(gca,'xticklabel',{[]}); 
subplot(6,13,[52+8 52+9 65+8 65+9]); gscatter([WMt{:,1};WMt{:,2}].*1000,[WMt{:,3};WMt{:,4}].*1000,gp);
grid on; axis('square'); xlabel('T1w'); legend('off'); axis([300 650 300 650])
subplot(6,13,[52+10 65+10]); plot(fliplr(K5{3}),'b','LineWidth',2); hold on; 
grid on; plot(fliplr(K5{4}),'r','LineWidth',2); set(gca,'xticklabel',{[]}); camroll(-90)

subplot(6,13,[39+11 39+12]); plot(K6{1},'b','LineWidth',2); hold on; plot(K6{2},'r','LineWidth',2); 
grid on; set(gca,'xticklabel',{[]}); 
subplot(6,13,[52+11 52+12 65+11 65+12]); gscatter([CSFt{:,1};CSFt{:,2}].*1000,[CSFt{:,3};CSFt{:,4}].*1000,gp);
grid on; axis('square'); xlabel('T1w'); legend('off'); axis([100 450 100 450])
subplot(6,13,[52+13 65+13]); plot(fliplr(K6{3}),'b','LineWidth',2); hold on; 
grid on; plot(fliplr(K6{4}),'r','LineWidth',2); set(gca,'xticklabel',{[]}); camroll(-90)

% table 1  
summary = table([TIVd_CI(1,1) TIVd_est(1) TIVd_CI(2,1); TIVt_CI(1,1) TIVt_est(1) TIVt_CI(2,1)],...
    [TIVd_CI(1,2) TIVd_est(2) TIVd_CI(2,2); TIVt_CI(1,2) TIVt_est(2) TIVt_CI(2,2)],...
    [TIVd_CI(1,3) TIVd_est(3) TIVd_CI(2,3); TIVt_CI(1,3) TIVt_est(3) TIVt_CI(2,3)],...
    [TIVd_CI(1,4) TIVd_est(4) TIVd_CI(2,4); TIVt_CI(1,4) TIVt_est(4) TIVt_CI(2,4)],...
    'RowNames',{'TIV discovery','TIV test'},'VariableNames',{'T1_nG1','T1_nG2','T12_nG1','T12_nG2'});
disp(summary); 

summary = table([CId_GM(1,1) GMd_est(1) CId_GM(2,1); CIt_GM(1,1) GMt_est(1) CIt_GM(2,1)],...
    [CId_GM(1,2) GMd_est(2) CId_GM(2,2); CIt_GM(1,2) GMt_est(2) CIt_GM(2,2)],...
    [CId_GM(1,3) GMd_est(3) CId_GM(2,3); CIt_GM(1,3) GMt_est(3) CIt_GM(2,3)],...
    [CId_GM(1,4) GMd_est(4) CId_GM(2,4); CIt_GM(1,4) GMt_est(4) CIt_GM(2,4)],...
    'RowNames',{'GM discovery','GM test'},'VariableNames',{'T1_nG1','T1_nG2','T12_nG1','T12_nG2'});
disp(summary); 

summary = table([CId_WM(1,1) WMd_est(1) CId_WM(2,1); CIt_WM(1,1) WMt_est(1) CIt_WM(2,1)],...
    [CId_WM(1,2) WMd_est(2) CId_WM(2,2); CIt_WM(1,2) WMt_est(2) CIt_WM(2,2)],...
    [CId_WM(1,3) WMd_est(3) CId_WM(2,3); CIt_WM(1,3) WMt_est(3) CIt_WM(2,3)],...
    [CId_WM(1,4) WMd_est(4) CId_WM(2,4); CIt_WM(1,4) WMt_est(4) CIt_WM(2,4)],...
    'RowNames',{'WM discovery','WM test'},'VariableNames',{'T1_nG1','T1_nG2','T12_nG1','T12_nG2'});
disp(summary); 

summary = table([CId_CSF(1,1) CSFd_est(1) CId_CSF(2,1); CIt_CSF(1,1) CSFt_est(1) CIt_CSF(2,1)],...
    [CId_CSF(1,2) CSFd_est(2) CId_CSF(2,2); CIt_CSF(1,2) CSFt_est(2) CIt_CSF(2,2)],...
    [CId_CSF(1,3) CSFd_est(3) CId_CSF(2,3); CIt_CSF(1,3) CSFt_est(3) CIt_CSF(2,3)],...
    [CId_CSF(1,4) CSFd_est(4) CId_CSF(2,4); CIt_CSF(1,4) CSFt_est(4) CIt_CSF(2,4)],...
    'RowNames',{'CSF discovery','CSF test'},'VariableNames',{'T1_nG1','T1_nG2','T12_nG1','T12_nG2'});
disp(summary); 

% -------------------------------------------------------------------------------
%% what is the total intracranial volume (TIV) for the four types of segmentation
% -------------------------------------------------------------------------------

% in the discovery set test main effects and interaction using a Hotelling
% test (repeated measure ANOVA) and multiple pair differences (alphav is adjusted
% using Hochberg step-up procedure)

result = rst_rep_anova_T2(TIVd,[],[2 2],1000,{'modality','n_gaussians'});
disp('-----');
warning('significant effect of the nb of gaussians %g ml and of modality %g ml, with no interaction',...
    mean(rst_trimmean(TIVd(:,[2 4])-TIVd(:,[1 3]))),mean(rst_trimmean(TIVd(:,[3 4])-TIVd(:,[1 2]))))
disp(result)
disp('-----')

% replication set - test for the same differences found as above using Bonferonni correction
Data1 = [TIVd(:,2)-TIVd(:,1), TIVd(:,4)-TIVd(:,3),...
    TIVd(:,3)-TIVd(:,1),TIVd(:,4)-TIVd(:,2)];
Data2 = [TIVt(:,2)-TIVt(:,1), TIVt(:,4)-TIVt(:,3),...
    TIVt(:,3)-TIVt(:,1),TIVt(:,4)-TIVt(:,2)];
[h,~,p] = rst_1ttest([mean(Data2(:,[1 2]),2),mean(Data2(:,[3 4]),2)],'estimator','trimmed mean','figure','off');
disp('-----');
warning('validation set confirms differences observed in the discovery set');
warning('adding 1 Gaussian increases TIV by %g',mean(rst_trimmean(TIVt(:,[2 4])-TIVt(:,[1 3]))));
warning('adding a T2 image decreases TIV by %g',mean(rst_trimmean(TIVt(:,[3 4])-TIVt(:,[1 2]))));

% estimates of differences for effect sizes
figure('Name','TIV differences'); subplot(1,2,1);
[TIVd_diff, CId_diff]   = rst_data_plot(Data1, 'estimator','trimmed mean','newfig','sub');
title('TIV differences discovery set','Fontsize',12); ylabel('TIV'); axis([0.1 5.9 -140 40]); subplot(1,2,2);
[TIVt_diff, CIt_diff]   = rst_data_plot(Data2, 'estimator','trimmed mean','newfig','sub');
title('TIV differences validation set','Fontsize',12); ylabel('TIV'); axis([0.1 5.9 -140 40]);

summary = table([CId_diff(1,1) TIVd_diff(1) CId_diff(2,1); CIt_diff(1,1) TIVt_diff(1) CIt_diff(2,1)],...
    [CId_diff(1,2) TIVd_diff(2) CId_diff(2,2); CIt_diff(1,2) TIVt_diff(2) CIt_diff(2,2)],...
    [CId_diff(1,3) TIVd_diff(3) CId_diff(2,3); CIt_diff(1,3) TIVt_diff(3) CIt_diff(2,3)],...
    [CId_diff(1,4) TIVd_diff(4) CId_diff(2,4); CIt_diff(1,4) TIVt_diff(4) CIt_diff(2,4)],...
    'RowNames',{'diff discovery','diff test'},'VariableNames',{'T1 G2-G1','T12 G2-G1','G1 T12-T1','G2 T12-T1'});
disp(summary); 

%% Where are the missing volumes?
% --------------------------------

SoftTissued = readtable(['nrudataset' filesep 'SoftTissue_volumes.csv'],'ReadRowNames',false);           
Skulld      = readtable(['nrudataset' filesep 'Skull_volumes.csv'],'ReadRowNames',false);           
Otherd      = readtable(['nrudataset' filesep 'Others_volumes.csv'],'ReadRowNames',false);   
SoftTissuet = readtable(['ds003653' filesep 'SoftTissue_volumes.csv'],'ReadRowNames',false);           
Skullt      = readtable(['ds003653' filesep 'Skull_volumes.csv'],'ReadRowNames',false);           
Othert      = readtable(['ds003653' filesep 'Others_volumes.csv'],'ReadRowNames',false);   

fprintf('Adding T2w leads to %g ml\n', (mean(TIVd(:,3)-TIVd(:,1)+TIVd(:,4)-TIVd(:,2))+...
    mean(TIVt(:,3)-TIVt(:,1)+TIVt(:,4)-TIVt(:,2)))/4)
fprintf('this is compensated by a change of %g ml in soft tissue\n', ...
    (mean(SoftTissued{:,3}-SoftTissued{:,1}+SoftTissued{:,4}-SoftTissued{:,2})+...
    mean(SoftTissuet{:,3}-SoftTissuet{:,1}+SoftTissuet{:,4}-SoftTissuet{:,2}))*250)
fprintf('this is compensated by a change of %g ml in bone tisue\n', ...
    (mean(Skulld{:,3}-Skulld{:,1}+Skulld{:,4}-Skulld{:,2})+...
    mean(Skullt{:,3}-Skullt{:,1}+Skullt{:,4}-Skullt{:,2}))*250)

fprintf('Adding a Gaussian leads to %g ml\n', (mean(TIVd(:,2)-TIVd(:,1)+TIVd(:,4)-TIVd(:,3))+...
    mean(TIVt(:,2)-TIVt(:,1)+TIVt(:,4)-TIVt(:,3)))/4)
fprintf('this is compensated by a change of %g ml in soft tissue\n', ...
    (mean(SoftTissued{:,2}-SoftTissued{:,1}+SoftTissued{:,4}-SoftTissued{:,4})+...
    mean(SoftTissuet{:,2}-SoftTissuet{:,1}+SoftTissuet{:,4}-SoftTissuet{:,4}))*250)
fprintf('this is compensated by a change of %g ml in bone tisue\n', ...
    (mean(Skulld{:,2}-Skulld{:,1}+Skulld{:,4}-Skulld{:,4})+...
    mean(Skullt{:,2}-Skullt{:,1}+Skullt{:,4}-Skullt{:,4}))*250)


% -------------------------------------------------------------
%% what are the volume differences for each brain tissue class 
% -------------------------------------------------------------

figure('Name','Tissue volume differences'); 

subplot(2,3,1);
Data = [GMd{:,2}-GMd{:,1}, GMd{:,4}-GMd{:,3}, GMd{:,3}-GMd{:,1}, GMd{:,4}-GMd{:,2}].*1000;
[GMd_diff, CIGMd_diff]   = rst_data_plot(Data, 'estimator','trimmed mean','newfig','sub');
title('GM differences discovery set','Fontsize',12);
subplot(2,3,4);
Data = [GMt{:,2}-GMt{:,1}, GMt{:,4}-GMt{:,3}, GMt{:,3}-GMt{:,1},GMt{:,4}-GMt{:,2}].*1000;
[GMt_diff, CIGMt_diff]   = rst_data_plot(Data, 'estimator','trimmed mean','newfig','sub');
title('GM differences validation set','Fontsize',12);

result = rst_rep_anova_T2(GMd{:,:},[],[2 2],1000,{'modality','n_gaussians'});
warning('significant effect of modality, nb of gaussians AND interaction') 
[GMdmeans, GMdCI] = rst_rep_anova_plot(GMd{:,:},ones(259,1),[2 2],3);
disp(result)
disp('-----')
[~,~,GMd_p,~,h1] = rst_multicompare(GMd{:,:}.*1000,[3 1; 4 2], 'estimator', 'trimmed mean','newfig','no');
fprintf('GM volumes differences by adding T2: %g for 1 Gaussians p=%g, %g for 2 Gaussians p=%g\n',...
    GMd_diff(3),GMd_p(1),GMd_diff(4),GMd_p(2))
[h1,CIx,GMd_p] = rst_1ttest((GMt{:,3}-GMt{:,1})-(GMt{:,4}-GMt{:,2})*1000,'trimmean');
fprintf('interation effect does not replicate %g\n',GMd_p)

figure(findobj( 'Type', 'Figure', 'Name', 'Tissue volume differences' ));
subplot(2,3,2);
Data = [WMd{:,2}-WMd{:,1}, WMd{:,4}-WMd{:,3}, WMd{:,3}-WMd{:,1}, WMd{:,4}-WMd{:,2}].*1000;
[WMd_diff, CIWMd_diff]   = rst_data_plot(Data, 'estimator','trimmed mean','newfig','sub');
title('WM differences discovery set','Fontsize',12);
subplot(2,3,5);
Data = [WMt{:,2}-WMt{:,1}, WMt{:,4}-WMt{:,3},  WMt{:,3}-WMt{:,1},WMt{:,4}-WMt{:,2}].*1000;
[WMt_diff, CIWMt_diff]   = rst_data_plot(Data, 'estimator','trimmed mean','newfig','sub');
title('WM differences validation set','Fontsize',12);

result = rst_rep_anova_T2(WMd{:,:},[],[2 2],1000,{'modality','n_gaussians'});
warning('significant effect of modality, nb of gaussians AND interaction') 
[WMdmeans, WMdCI] = rst_rep_anova_plot(WMd{:,:},ones(259,1),[2 2],3);
disp(result)
disp('-----')
[~,~,WMd_p,~,h2] = rst_multicompare(WMd{:,:}.*1000,[3 1; 4 2], 'estimator', 'trimmed mean','newfig','no');
fprintf('WM volumes differences by adding T2: %g for 1 Gaussians p=%g, %g for 2 Gaussians p=%g\n', ...
    WMd_diff(3),WMd_p(1),WMd_diff(4),WMd_p(2))
[h2,CIx,WMd_p] = rst_1ttest((WMt{:,3}-WMt{:,1})-(WMt{:,4}-WMt{:,2})*1000,'trimmean');
fprintf('interation effect does not replicate %g\n',WMd_p)

figure(findobj( 'Type', 'Figure', 'Name', 'Tissue volume differences' ));
subplot(2,3,3);
Data = [CSFd{:,2}-CSFd{:,1}, CSFd{:,4}-CSFd{:,3}, CSFd{:,3}-CSFd{:,1}, CSFd{:,4}-CSFd{:,2}].*1000;
[CSFd_diff, CICSFd_diff]   = rst_data_plot(Data, 'estimator','trimmed mean','newfig','sub');
title('CSF differences discovery set','Fontsize',12);
subplot(2,3,6);
Data = [CSFt{:,2}-CSFt{:,1}, CSFt{:,4}-CSFt{:,3}, CSFt{:,3}-CSFt{:,1},CSFt{:,4}-CSFt{:,2}].*1000;
[CSFt_diff, CICSFt_diff]   = rst_data_plot(Data, 'estimator','trimmed mean','newfig','sub');
title('CSF differences validation set','Fontsize',12);

result = rst_rep_anova_T2(CSFd{:,:},[],[2 2],1000,{'modality','n_gaussians'});
warning('significant effect of modality, nb of gaussians AND interaction') 
[CSFdmeans, CSFdCI] = rst_rep_anova_plot(CSFd{:,:},ones(259,1),[2 2],3);
disp(result)
disp('-----')
[~,~,CSFd_p,~,h3] = rst_multicompare(CSFd{:,:}.*1000,[3 1; 4 2], 'estimator', 'trimmed mean','newfig','no');
fprintf('CSF volumes differences by adding T2: %g for 1 Gaussians p=%g, %g for 2 Gaussians p=%g\n', ...
    CSFd_diff(3),CSFd_p(1),CSFd_diff(4),CSFd_p(2))
[h3,CIx,CSFd_p] = rst_1ttest((WMt{:,3}-WMt{:,1})-(WMt{:,4}-WMt{:,2})*1000,'trimmean');
fprintf('interation effect does not replicate %g\n',CSFd_p)

% summary tables
% --------------
summary = table([CIGMd_diff(1,1) GMd_diff(1) CIGMd_diff(2,1); CIGMt_diff(1,1) GMt_diff(1) CIGMt_diff(2,1)],...
    [CIGMd_diff(1,2) GMd_diff(2) CIGMd_diff(2,2); CIGMt_diff(1,2) GMt_diff(2) CIGMt_diff(2,2)],...
    [CIGMd_diff(1,3) GMd_diff(3) CIGMd_diff(2,3); CIGMt_diff(1,3) GMt_diff(3) CIGMt_diff(2,3)],...
    [CIGMd_diff(1,4) GMd_diff(4) CIGMd_diff(2,4); CIGMt_diff(1,4) GMt_diff(4) CIGMt_diff(2,4)],...
    'RowNames',{'diff GM discovery','diff GM test'},'VariableNames',{'T1 G2-G1','T12 G2-G1','G1 T12-T1','G2 T12-T1'});
disp(summary); 

summary = table([CIWMd_diff(1,1) WMd_diff(1) CIWMd_diff(2,1); CIWMt_diff(1,1) WMt_diff(1) CIWMt_diff(2,1)],...
    [CIWMd_diff(1,2) WMd_diff(2) CIWMd_diff(2,2); CIWMt_diff(1,2) WMt_diff(2) CIWMt_diff(2,2)],...
    [CIWMd_diff(1,3) WMd_diff(3) CIWMd_diff(2,3); CIWMt_diff(1,3) WMt_diff(3) CIWMt_diff(2,3)],...
    [CIWMd_diff(1,4) WMd_diff(4) CIWMd_diff(2,4); CIWMt_diff(1,4) WMt_diff(4) CIWMt_diff(2,4)],...
    'RowNames',{'diff WM discovery','diff WM test'},'VariableNames',{'T1 G2-G1','T12 G2-G1','G1 T12-T1','G2 T12-T1'});disp(summary); 

summary = table([CICSFd_diff(1,1) CSFd_diff(1) CICSFd_diff(2,1); CICSFt_diff(1,1) CSFt_diff(1) CICSFt_diff(2,1)],...
    [CICSFd_diff(1,2) CSFd_diff(2) CICSFd_diff(2,2); CICSFt_diff(1,2) CSFt_diff(2) CICSFt_diff(2,2)],...
    [CICSFd_diff(1,3) CSFd_diff(3) CICSFd_diff(2,3); CICSFt_diff(1,3) CSFt_diff(3) CICSFt_diff(2,3)],...
    [CICSFd_diff(1,4) CSFd_diff(4) CICSFd_diff(2,4); CICSFt_diff(1,4) CSFt_diff(4) CICSFt_diff(2,4)],...
    'RowNames',{'diff CSF discovery','diff CSF test'},'VariableNames',{'T1 G2-G1','T12 G2-G1','G1 T12-T1','G2 T12-T1'});
disp(summary); 

% tissue volume changes can be seen in classifiying non-brain tissue
% located inside the brain too (not just solft tissue and bones)

% For vessels, we computed the % of voxels being not Grey (<0.1) or Grey
% (>0.9), not white (<0.1) or white (>0.9), and not csf (<0.1) or csf (>0.9). 
% We summarize this here by using the ratio, if the ratio is bigger than 1
% it indicates more voxels seen as not from that tissue than from that tissue
% % and conversely, the higher that ratio the better.

GMd       = readtable(['nrudataset' filesep 'GreyMatter_vessels.csv'],'ReadRowNames',false);  % High probability of Grey matter in vessels
WMd       = readtable(['nrudataset' filesep 'WhiteMatter_vessels.csv'],'ReadRowNames',false);           
CSFd      = readtable(['nrudataset' filesep 'CSF_vessels.csv'],'ReadRowNames',false);   
notGMd    = readtable(['nrudataset' filesep 'Not_GreyMatter_vessels.csv'],'ReadRowNames',false);  % Low probability of Grey matter in vessels
notWMd    = readtable(['nrudataset' filesep 'Not_WhiteMatter_vessels.csv'],'ReadRowNames',false);           
notCSFd   = readtable(['nrudataset' filesep 'Not_CSF_vessels.csv'],'ReadRowNames',false);   
GMratiod  = [mean(notGMd{:,[1 2]} ./ GMd{:,[1 2]},2) mean(notGMd{:,[3 4]} ./ GMd{:,[3 4]},2)];
WMratiod  = [mean(notWMd{:,[1 2]} ./ WMd{:,[1 2]},2) mean(notWMd{:,[3 4]} ./ WMd{:,[3 4]},2)];
CSFratiod = [mean(notCSFd{:,[1 2]}./CSFd{:,[1 2]},2) mean(notCSFd{:,[3 4]}./CSFd{:,[3 4]},2)];

GMt       = readtable(['ds003653' filesep 'GreyMatter_vessels.csv'],'ReadRowNames',false);           
WMt       = readtable(['ds003653' filesep 'WhiteMatter_vessels.csv'],'ReadRowNames',false);           
CSFt      = readtable(['ds003653' filesep 'CSF_vessels.csv'],'ReadRowNames',false);
notGMt    = readtable(['ds003653' filesep 'Not_GreyMatter_vessels.csv'],'ReadRowNames',false);           
notWMt    = readtable(['ds003653' filesep 'Not_WhiteMatter_vessels.csv'],'ReadRowNames',false);           
notCSFt   = readtable(['ds003653' filesep 'Not_CSF_vessels.csv'],'ReadRowNames',false);
GMratiot  = [mean(notGMt{:,[1 2]} ./ GMt{:,[1 2]},2) mean(notGMt{:,[3 4]} ./ GMt{:,[3 4]},2)];
WMratiot  = [mean(notWMt{:,[1 2]} ./ WMt{:,[1 2]},2) mean(notWMt{:,[3 4]} ./ WMt{:,[3 4]},2)];
CSFratiot = [mean(notCSFt{:,[1 2]}./CSFt{:,[1 2]},2) mean(notCSFt{:,[3 4]}./CSFt{:,[3 4]},2)];

% Discovery set
figure('Name','Vessels tissue attribution'); 
subplot(3,4,1); boxchart([mean(GMd{:,[1 2]},2) mean(GMd{:,[3 4]},2)])
hold on; plot(median([mean(GMd{:,[1 2]},2) mean(GMd{:,[3 4]},2)]),'.','MarkerSize',25)
xticklabels({'T1','T1 and T2'}); 
title('Grey Matter discovery set','Fontsize',12); 
ylabel('% of voxels as GM'); grid on; box on; ylim([10 40])
subplot(3,4,5); boxchart([mean(WMd{:,[1 2]},2) mean(WMd{:,[3 4]},2)])
hold on; plot(median([mean(WMd{:,[1 2]},2) mean(WMd{:,[3 4]},2)]),'.','MarkerSize',25)
xticklabels({'T1','T1 and T2'}); 
title('White Matter discovery set','Fontsize',12); 
ylabel('% of voxels as WM'); grid on; box on; ylim([0 4])
subplot(3,4,9); boxchart([mean(CSFd{:,[1 2]},2) mean(CSFd{:,[3 4]},2)])
hold on; plot(median([mean(CSFd{:,[1 2]},2) mean(CSFd{:,[3 4]},2)]),'.','MarkerSize',25)
xticklabels({'T1','T1 and T2'}); 
title('CSF discovery set','Fontsize',12);
ylabel('% of voxels as CSF'); grid on; box on; ylim([0 30])

% replication set
subplot(3,4,2); boxchart([mean(GMt{:,[1 2]},2) mean(GMt{:,[3 4]},2)])
hold on; plot(median([mean(GMt{:,[1 2]},2) mean(GMt{:,[3 4]},2)]),'.','MarkerSize',25)
xticklabels({'T1','T1 and T2'}); 
title('Grey Matter validation set','Fontsize',12); 
ylabel('% of voxels as GM'); grid on; box on; ylim([10 40])
subplot(3,4,6); boxchart([mean(WMt{:,[1 2]},2) mean(WMt{:,[3 4]},2)])
hold on; plot(median([mean(WMt{:,[1 2]},2) mean(WMt{:,[3 4]},2)]),'.','MarkerSize',25)
xticklabels({'T1','T1 and T2'}); 
title('White Matter validation set','Fontsize',12); 
ylabel('% of voxels as WM'); grid on; box on; ylim([0 4])
subplot(3,4,10); boxchart([mean(CSFt{:,[1 2]},2) mean(CSFt{:,[3 4]},2)])
hold on; plot(median([mean(CSFt{:,[1 2]},2) mean(CSFt{:,[3 4]},2)]),'.','MarkerSize',25)
xticklabels({'T1','T1 and T2'}); 
title('CSF validation set','Fontsize',12);
ylabel('% of voxels as CSF'); grid on; box on; ylim([0 30])

% Discovery set
subplot(3,4,3); boxchart([mean(notGMd{:,[1 2]},2) mean(notGMd{:,[3 4]},2)])
hold on; plot(median([mean(notGMd{:,[1 2]},2) mean(notGMd{:,[3 4]},2)]),'.','MarkerSize',25)
xticklabels({'T1','T1 and T2'}); 
title('Not Grey Matter discovery set','Fontsize',12); 
ylabel('% of voxels as not GM'); grid on; box on; ylim([20 60])
subplot(3,4,7); boxchart([mean(notWMd{:,[1 2]},2) mean(notWMd{:,[3 4]},2)])
hold on; plot(median([mean(notWMd{:,[1 2]},2) mean(notWMd{:,[3 4]},2)]),'.','MarkerSize',25)
xticklabels({'T1','T1 and T2'}); 
title('Not White Matter discovery set','Fontsize',12); 
ylabel('% of voxels as not WM'); grid on; box on; ylim([80 95])
subplot(3,4,11); boxchart([mean(notCSFd{:,[1 2]},2) mean(notCSFd{:,[3 4]},2)])
hold on; plot(median([mean(notCSFd{:,[1 2]},2) mean(notCSFd{:,[3 4]},2)]),'.','MarkerSize',25)
xticklabels({'T1','T1 and T2'}); 
title('Not CSF discovery set','Fontsize',12);
ylabel('% of voxels as not CSF'); grid on; box on; ylim([30 70])

% validation set
subplot(3,4,4); boxchart([mean(notGMt{:,[1 2]},2) mean(notGMt{:,[3 4]},2)])
hold on; plot(median([mean(notGMt{:,[1 2]},2) mean(notGMt{:,[3 4]},2)]),'.','MarkerSize',25)
xticklabels({'T1','T1 and T2'}); 
title('Not Grey Matter validation set','Fontsize',12); 
ylabel('% of voxels as not GM'); grid on; box on; ylim([20 60])
subplot(3,4,8); boxchart([mean(notWMt{:,[1 2]},2) mean(notWMt{:,[3 4]},2)])
hold on; plot(median([mean(notWMt{:,[1 2]},2) mean(notWMt{:,[3 4]},2)]),'.','MarkerSize',25)
xticklabels({'T1','T1 and T2'}); 
title('Not White Matter validation set','Fontsize',12); 
ylabel('% of voxels as not WM'); grid on; box on; ylim([80 95])
subplot(3,4,12); boxchart([mean(notCSFt{:,[1 2]},2) mean(notCSFt{:,[3 4]},2)])
hold on; plot(median([mean(notCSFt{:,[1 2]},2) mean(notCSFt{:,[3 4]},2)]),'.','MarkerSize',25)
xticklabels({'T1','T1 and T2'}); 
title('Not CSF validation set','Fontsize',12);
ylabel('% of voxels as not CSF'); grid on; box on; ylim([30 70])

% check that summary stats makes sense
GMdmeans = mean(GMd{:,:})+mean(notGMd{:,:});
WMdmeans = mean(WMd{:,:})+mean(notWMd{:,:});
CSFdmeans = mean(CSFd{:,:})+mean(notCSFd{:,:});
GMtmeans = mean(GMt{:,:})+mean(notGMt{:,:});
WMtmeans = mean(WMt{:,:})+mean(notWMt{:,:});
CSFtmeans = mean(CSFt{:,:})+mean(notCSFt{:,:});
warning('While not capturing the full range of probabilities, threshoding tissues at 0.1 and 0.9, captured more than half of all vessel voxels')
warning('Over all voxels containing vessels, %g%% [min %g max %g] for discovery data and %g%% [min %g max %g] for validation data were in the lower and higher GM centiles', ...
    mean(GMdmeans),min(GMdmeans),max(GMdmeans),mean(GMtmeans),min(GMtmeans),max(GMtmeans))
warning('%g%% [min %g max %g]  and %g%% [min %g max %g] were in the lower and higher WM centiles',...
    mean(WMdmeans),min(WMdmeans),max(WMdmeans),mean(WMtmeans),min(WMtmeans),max(WMtmeans))
warning('and %g%% [min %g max %g]  and %g%% [min %g max %g] were in the lower and higher CSF centiles',...
    mean(CSFdmeans),min(CSFdmeans),max(CSFdmeans),mean(CSFtmeans),min(CSFtmeans),max(CSFtmeans))
warning('Simply looking at the percentages, results indicates that voxels containing vessels are moslty seen as a mixture of GM and CSF')

% test for differences in ratios
figure('Name','ratio tests'); subplot(1,2,1)
data = [diff(GMratiod,1,2) diff(WMratiod,1,2) diff(CSFratiod,1,2)];
[hd,CId,pd] = rst_1ttest(data,'trimmean','newfig','no'); title('Discovery set ratios')
warning('In the discovery set, adding the T2w image leads to increase the number of voxels not seen as belonging to GM. WM and CSF (p<.0017)')
data = [diff(GMratiot,1,2) diff(WMratiot,1,2) diff(CSFratiot,1,2)]; subplot(1,2,2)
[ht,CIt,pt] = rst_1ttest(data,'trimmean','newfig','no'); title('Validation set ratios')
warning('While in the validation set, adding the T2w image also leads to increase the number of voxels not seen as belonging to GM. WM and CSF (p<.0017)')
warning('the larger change is seen for CSF while in the discovery set this was for WM')



%% Multivariate analysis of volumes
% ---------------------------------

% let's simply count how may people change in each direction
% GM+WM+CSF+ // GM+WM+CSF- // GM+WM-CSF+ // GM+WM-CSF- 
% GM-WM+CSF+ // GM-WM+CSF- // GM-WM-CSF+ // GM-WM-CSF- 

gp_11 = (GMd{:,3}-GMd{:,1}>0).*(WMd{:,3}-WMd{:,1}>0).*(CSFd{:,3}-CSFd{:,1}>0);
gp_21 = (GMd{:,3}-GMd{:,1}>0).*(WMd{:,3}-WMd{:,1}>0).*(CSFd{:,3}-CSFd{:,1}<0);
gp_31 = (GMd{:,3}-GMd{:,1}>0).*(WMd{:,3}-WMd{:,1}<0).*(CSFd{:,3}-CSFd{:,1}>0);
gp_41 = (GMd{:,3}-GMd{:,1}>0).*(WMd{:,3}-WMd{:,1}<0).*(CSFd{:,3}-CSFd{:,1}<0);
gp_51 = (GMd{:,3}-GMd{:,1}<0).*(WMd{:,3}-WMd{:,1}>0).*(CSFd{:,3}-CSFd{:,1}>0);
gp_61 = (GMd{:,3}-GMd{:,1}<0).*(WMd{:,3}-WMd{:,1}>0).*(CSFd{:,3}-CSFd{:,1}<0);
gp_71 = (GMd{:,3}-GMd{:,1}<0).*(WMd{:,3}-WMd{:,1}<0).*(CSFd{:,3}-CSFd{:,1}>0);
gp_81 = (GMd{:,3}-GMd{:,1}<0).*(WMd{:,3}-WMd{:,1}<0).*(CSFd{:,3}-CSFd{:,1}<0);

gp_12 = (GMd{:,4}-GMd{:,2}>0).*(WMd{:,4}-WMd{:,2}>0).*(CSFd{:,4}-CSFd{:,2}>0);
gp_22 = (GMd{:,4}-GMd{:,2}>0).*(WMd{:,4}-WMd{:,2}>0).*(CSFd{:,4}-CSFd{:,2}<0);
gp_32 = (GMd{:,4}-GMd{:,2}>0).*(WMd{:,4}-WMd{:,2}<0).*(CSFd{:,4}-CSFd{:,2}>0);
gp_42 = (GMd{:,4}-GMd{:,2}>0).*(WMd{:,4}-WMd{:,2}<0).*(CSFd{:,4}-CSFd{:,2}<0);
gp_52 = (GMd{:,4}-GMd{:,2}<0).*(WMd{:,4}-WMd{:,2}>0).*(CSFd{:,4}-CSFd{:,2}>0);
gp_62 = (GMd{:,4}-GMd{:,2}<0).*(WMd{:,4}-WMd{:,2}>0).*(CSFd{:,4}-CSFd{:,2}<0);
gp_72 = (GMd{:,4}-GMd{:,2}<0).*(WMd{:,4}-WMd{:,2}<0).*(CSFd{:,4}-CSFd{:,2}>0);
gp_82 = (GMd{:,4}-GMd{:,2}<0).*(WMd{:,4}-WMd{:,2}<0).*(CSFd{:,4}-CSFd{:,2}<0);

summary = table(mean([gp_11 gp_21 gp_22 gp_22]',2)*100, mean([gp_31 gp_41 gp_32 gp_42]',2)*100', ...
    mean([gp_51 gp_61 gp_52 gp_62]',2)*100, mean([gp_71 gp_81 gp_72 gp_82]',2)*100, ...
    'RowNames',{'CSF+ 1 Gaussian','CSF- 1 Gaussian', 'CSF+ 2 Gaussians', 'CSF- 2 Gaussians'},...
    'VariableNames',{'GM+WM+','GM+WM-','GM-WM+','GM-WM-'});
warning('Discovery dataset'); disp(summary); 

gp2_11 = (GMt{:,3}-GMt{:,1}>0).*(WMt{:,3}-WMt{:,1}>0).*(CSFt{:,3}-CSFt{:,1}>0);
gp2_21 = (GMt{:,3}-GMt{:,1}>0).*(WMt{:,3}-WMt{:,1}>0).*(CSFt{:,3}-CSFt{:,1}<0);
gp2_31 = (GMt{:,3}-GMt{:,1}>0).*(WMt{:,3}-WMt{:,1}<0).*(CSFt{:,3}-CSFt{:,1}>0);
gp2_41 = (GMt{:,3}-GMt{:,1}>0).*(WMt{:,3}-WMt{:,1}<0).*(CSFt{:,3}-CSFt{:,1}<0);
gp2_51 = (GMt{:,3}-GMt{:,1}<0).*(WMt{:,3}-WMt{:,1}>0).*(CSFt{:,3}-CSFt{:,1}>0);
gp2_61 = (GMt{:,3}-GMt{:,1}<0).*(WMt{:,3}-WMt{:,1}>0).*(CSFt{:,3}-CSFt{:,1}<0);
gp2_71 = (GMt{:,3}-GMt{:,1}<0).*(WMt{:,3}-WMt{:,1}<0).*(CSFt{:,3}-CSFt{:,1}>0);
gp2_81 = (GMt{:,3}-GMt{:,1}<0).*(WMt{:,3}-WMt{:,1}<0).*(CSFt{:,3}-CSFt{:,1}<0);

gp2_12 = (GMt{:,4}-GMt{:,2}>0).*(WMt{:,4}-WMt{:,2}>0).*(CSFt{:,4}-CSFt{:,2}>0);
gp2_22 = (GMt{:,4}-GMt{:,2}>0).*(WMt{:,4}-WMt{:,2}>0).*(CSFt{:,4}-CSFt{:,2}<0);
gp2_32 = (GMt{:,4}-GMt{:,2}>0).*(WMt{:,4}-WMt{:,2}<0).*(CSFt{:,4}-CSFt{:,2}>0);
gp2_42 = (GMt{:,4}-GMt{:,2}>0).*(WMt{:,4}-WMt{:,2}<0).*(CSFt{:,4}-CSFt{:,2}<0);
gp2_52 = (GMt{:,4}-GMt{:,2}<0).*(WMt{:,4}-WMt{:,2}>0).*(CSFt{:,4}-CSFt{:,2}>0);
gp2_62 = (GMt{:,4}-GMt{:,2}<0).*(WMt{:,4}-WMt{:,2}>0).*(CSFt{:,4}-CSFt{:,2}<0);
gp2_72 = (GMt{:,4}-GMt{:,2}<0).*(WMt{:,4}-WMt{:,2}<0).*(CSFt{:,4}-CSFt{:,2}>0);
gp2_82 = (GMt{:,4}-GMt{:,2}<0).*(WMt{:,4}-WMt{:,2}<0).*(CSFt{:,4}-CSFt{:,2}<0);

summary = table(mean([gp2_11 gp2_21 gp2_22 gp2_22]',2)*100, mean([gp2_31 gp2_41 gp2_32 gp2_42]',2)*100', ...
    mean([gp2_51 gp2_61 gp2_52 gp2_62]',2)*100, mean([gp2_71 gp2_81 gp2_72 gp2_82]',2)*100, ...
    'RowNames',{'CSF+ 1 Gaussian','CSF- 1 Gaussian', 'CSF+ 2 Gaussians', 'CSF- 2 Gaussians'},...
    'VariableNames',{'GM+WM+','GM+WM-','GM-WM+','GM-WM-'});
warning('Validation dataset'); disp(summary); 

% Now do a hard clustering
% ------------------------

Datad = [GMd{:,3}-GMd{:,1},WMd{:,3}-WMd{:,1},CSFd{:,3}-CSFd{:,1}];
[BICS,BESTMODEL] = mbclust(Datad,8);
figure('Name','Guassian Mixture Modelling')
subplot(2,5,1); plotbic(BICS); grid on; box on; axis square; axis([0.5 8.5 3800 4900]);
[class1,uncertainty] = mixclass(Datad,BESTMODEL.pies,BESTMODEL.mus,BESTMODEL.vars);
C   = zeros(length(class1),3);
C(class1==1,1) = 1; % red
C(class1==2,2) = 1; % green
C(class1==3,3) = 1; % blue
subplot(2,5,2);
scatter3(Datad(:,1),Datad(:,3),Datad(:,2),30,C.*(1-uncertainty)','filled'); 
xlabel('GM'); ylabel('CSF'); zlabel('WM'); axis square; axis([-0.1 0.05 -0.15 0.1 -0.06 0.01])
title(sprintf('Discovery set, 3 clusters\n mean error: %g',mean(uncertainty)))

Datat = [GMt{:,3}-GMt{:,1},WMt{:,3}-WMt{:,1},CSFt{:,3}-CSFt{:,1}];
[class,uncertainty] = mixclass(Datat,BESTMODEL.pies,BESTMODEL.mus,BESTMODEL.vars);
C   = zeros(length(class),3);
C(class==1,1) = 1; % red
C(class==2,2) = 1; % green
C(class==3,3) = 1; % blue
subplot(2,5,3);
scatter3(Datat(:,1),Datat(:,3),Datat(:,2),30,C.*(1-uncertainty)','filled'); 
xlabel('GM'); ylabel('CSF'); zlabel('WM'); axis square; axis([-0.1 0.05 -0.15 0.1 -0.06 0.01])
title(sprintf('Validation set, same model\n mean error: %g',mean(uncertainty)))

[BICS,BESTMODEL] = mbclust(Datat,8);
subplot(2,5,4); plotbic(BICS); grid on; box on; axis square; axis([0.5 8.5 1450 1800]);
[class,uncertainty] = mixclass(Datat,BESTMODEL.pies,BESTMODEL.mus,BESTMODEL.vars);
C = zeros(length(class),3);
colours = rst_colour_maps(8);
for c=1:8; C(class==c,:) = repmat(colours(c,:),sum(class==c),1); end
subplot(2,5,5);
scatter3(Datat(:,1),Datat(:,3),Datat(:,2),30,C.*(1-uncertainty)','filled'); 
xlabel('GM'); ylabel('CSF'); zlabel('WM'); axis square; axis([-0.1 0.05 -0.15 0.1 -0.06 0.01])
title(sprintf('Validation set, 8 clusters\n mean error: %g',mean(uncertainty)))

Datad = [GMd{:,4}-GMd{:,2},WMd{:,4}-WMd{:,2},CSFd{:,4}-CSFd{:,2}];
[BICS,BESTMODEL] = mbclust(Datad,8); clustering2Gd = BESTMODEL;
subplot(2,5,6); plotbic(BICS); grid on; box on; axis square; axis([0.5 8.5 3800 4900]);
[class2,uncertainty] = mixclass(Datad,BESTMODEL.pies,BESTMODEL.mus,BESTMODEL.vars);
subplot(2,5,7);
C   = zeros(length(class2),3);
C(class2==3,1) = 1; % red
C(class2==2,2) = 1; % green
C(class2==1,3) = 1; % blue
scatter3(Datad(:,1),Datad(:,3),Datad(:,2),30,C.*(1-uncertainty)','filled'); 
xlabel('GM'); ylabel('CSF'); zlabel('WM'); axis square; axis([-0.1 0.05 -0.15 0.1 -0.005 0.005])
title(sprintf('Discovery set, 3 clusters\n mean error: %g',mean(uncertainty)))

Datat = [GMt{:,4}-GMt{:,2},WMt{:,4}-WMt{:,2},CSFt{:,4}-CSFt{:,2}];
[class,uncertainty] = mixclass(Datat,BESTMODEL.pies,BESTMODEL.mus,BESTMODEL.vars);
subplot(2,5,8);
C   = zeros(length(class),3);
C(class==3,1) = 1; % red
C(class==2,2) = 1; % green
C(class==1,3) = 1; % blue
scatter3(Datat(:,1),Datat(:,3),Datat(:,2),30,C.*(1-uncertainty)','filled'); 
xlabel('GM'); ylabel('CSF'); zlabel('WM'); axis square; axis([-0.1 0.05 -0.15 0.1 -0.005 0.005])
title(sprintf('Validation set, same model\n mean error: %g',mean(uncertainty)))

[BICS,BESTMODEL,ALLMODELS] = mbclust(Datat,8); clustering2Gt = BESTMODEL;
subplot(2,5,9); plotbic(BICS); grid on; box on; axis square; axis([0.5 8.5 1450 1800]);
[class,uncertainty] = mixclass(Datat,BESTMODEL.pies,BESTMODEL.mus,BESTMODEL.vars);
subplot(2,5,10);
C   = zeros(length(class),3);
C(class==2,3) = 1; % blue
C(class==1,2) = 1; % green
scatter3(Datat(:,1),Datat(:,3),Datat(:,2),30,C.*(1-uncertainty)','filled'); 
xlabel('GM'); ylabel('CSF'); zlabel('WM'); axis square; axis([-0.1 0.05 -0.15 0.1 -0.005 0.005])
title(sprintf('Validation set, 2 clusters\n mean error: %g',mean(uncertainty)))

% we can see that some subjects are different - worth tracking them
% ie subjects 57, 116, 126, 144
outlier_class = intersect(find(class2==3),find(class1==1));
[~,~,up1] = rst_trimci((CSFd{:,3}-CSFd{:,1})*1000); 
[~,~,up2] = rst_trimci((CSFd{:,4}-CSFd{:,2})*1000);
fprintf('those 4 subjects have CSF changes of %g and %g\nwhile group estimates upper bound are %g and %g',...
    mean(CSFd{outlier_class,3}-CSFd{outlier_class,1})*1000, ...
    mean(CSFd{outlier_class,4}-CSFd{outlier_class,2})*1000, up1,up2)

% same main clusters
fprintf('Main clusters are similar when using 2 Gaussians mu=[%g %g %g] vs [%g %g %g]\n',...
    clustering2Gt.mus(:,1)*1000, clustering2Gd.mus(:,2)*1000)
fprintf('sigma= [%g %g %g] vs [%g %g %g]\n', diag(clustering2Gd.vars(:,:,2))*1000, ...
    diag(clustering2Gt.vars(:,:,1))*1000)


% -------------------------------------------------------------------------
%% What does this all mean in terms of similarlity/differences among tissues
% -------------------------------------------------------------------------

%% Dunn index
GMd  = readtable(['nrudataset' filesep 'GrayMatter_DunnIndexes.csv'],'ReadRowNames',false);           
WMd  = readtable(['nrudataset' filesep 'WhiteMatter_DunnIndexes.csv'],'ReadRowNames',false);           
CSFd = readtable(['nrudataset' filesep 'CSF_DunnIndexes.csv'],'ReadRowNames',false);           
GMt  = readtable(['ds003653' filesep 'GrayMatter_DunnIndexes.csv'],'ReadRowNames',false);           
WMt  = readtable(['ds003653' filesep 'WhiteMatter_DunnIndexes.csv'],'ReadRowNames',false);           
CSFt = readtable(['ds003653' filesep 'CSF_DunnIndexes.csv'],'ReadRowNames',false);           

[GMd_est, CId_GM]   = rst_trimmean(GMd{:,:});
[WMd_est, CId_WM]   = rst_trimmean(WMd{:,:});
[CSFd_est, CId_CSF] = rst_trimmean(CSFd{:,:});
 
TrimmedMeans = [GMd_est; WMd_est; CSFd_est];
LowerConfs   = [CId_GM(1,:); CId_WM(1,:); CId_CSF(1,:)];
HigherConfs  = [CId_GM(2,:); CId_WM(2,:); CId_CSF(2,:)];

T1_nG1  = [LowerConfs(:,1) TrimmedMeans(:,1) HigherConfs(:,1)];
T1_nG2  = [LowerConfs(:,2) TrimmedMeans(:,2) HigherConfs(:,2)];
T12_nG1 = [LowerConfs(:,3) TrimmedMeans(:,3) HigherConfs(:,3)];
T12_nG2 = [LowerConfs(:,4) TrimmedMeans(:,4) HigherConfs(:,4)];

DI_table = table(T1_nG1,T1_nG2,T12_nG1,T12_nG2,...
 'RowNames',{'GM','WM','CSF'});
disp(DI_table)

% validation set
[GMt_est, CIt_GM]   = rst_trimmean(GMt{:,:});
[WMt_est, CIt_WM]   = rst_trimmean(WMt{:,:});
[CSFt_est, CIt_CSF] = rst_trimmean(CSFt{:,:});
 
TrimmedMeans = [GMt_est; WMt_est; CSFt_est];
LowerConfs   = [CIt_GM(1,:); CIt_WM(1,:); CIt_CSF(1,:)];
HigherConfs  = [CIt_GM(2,:); CIt_WM(2,:); CIt_CSF(2,:)];

T1_nG1  = [LowerConfs(:,1) TrimmedMeans(:,1) HigherConfs(:,1)];
T1_nG2  = [LowerConfs(:,2) TrimmedMeans(:,2) HigherConfs(:,2)];
T12_nG1 = [LowerConfs(:,3) TrimmedMeans(:,3) HigherConfs(:,3)];
T12_nG2 = [LowerConfs(:,4) TrimmedMeans(:,4) HigherConfs(:,4)];

DI_table = table(T1_nG1,T1_nG2,T12_nG1,T12_nG2,...
 'RowNames',{'GM','WM','CSF'});
disp(DI_table)

disp("--------")
disp("There are not any evidence for")
disp("the conditions changing the Dunn Index")
disp("This tells that the tissues do not change value ranges")
disp("between the different conditions")
disp("--------")

%% Entropy
GMd  = readtable(['nrudataset' filesep 'GrayMatter_entropy.csv'],'ReadRowNames',false);           
WMd  = readtable(['nrudataset' filesep 'WhiteMatter_entropy.csv'],'ReadRowNames',false);           
CSFd = readtable(['nrudataset' filesep 'CSF_entropy.csv'],'ReadRowNames',false);           
GMt  = readtable(['ds003653' filesep 'GrayMatter_entropy.csv'],'ReadRowNames',false);           
WMt  = readtable(['ds003653' filesep 'WhiteMatter_entropy.csv'],'ReadRowNames',false);           
CSFt = readtable(['ds003653' filesep 'CSF_entropy.csv'],'ReadRowNames',false);   

figure('Name','Tissue Entropy'); 
subplot(3,4,1);
[GMd_est, CId_GM]   = rst_data_plot(GMd{:,:}, 'estimator','trimmed mean','newfig','sub');
title('Grey Matter discovery set','Fontsize',12); ylabel('GM Entropy'); subplot(3,4,5);
[WMd_est, CId_WM]   = rst_data_plot(WMd{:,:}, 'estimator','trimmed mean','newfig','sub');
title('White Matter discovery set','Fontsize',12); ylabel('WM Entropy'); subplot(3,4,9);
[CSFd_est, CId_CSF] = rst_data_plot(CSFd{:,:}, 'estimator','trimmed mean','newfig','sub');
title('CSF discovery set','Fontsize',12); ylabel('CSF Entropy');
 
TrimmedMeans = [GMd_est; WMd_est; CSFd_est];
LowerConfs   = [CId_GM(1,:); CId_WM(1,:); CId_CSF(1,:)];
HigherConfs  = [CId_GM(2,:); CId_WM(2,:); CId_CSF(2,:)];

T1_nG1  = [LowerConfs(:,1) TrimmedMeans(:,1) HigherConfs(:,1)];
T1_nG2  = [LowerConfs(:,2) TrimmedMeans(:,2) HigherConfs(:,2)];
T12_nG1 = [LowerConfs(:,3) TrimmedMeans(:,3) HigherConfs(:,3)];
T12_nG2 = [LowerConfs(:,4) TrimmedMeans(:,4) HigherConfs(:,4)];

DI_table = table(T1_nG1,T1_nG2,T12_nG1,T12_nG2,...
 'RowNames',{'GM','WM','CSF'});

subplot(3,4,2);
[GMd_diff,GMd_dCI,GMd_p,GMd_alphav,h1] = rst_multicompare(GMd{:,:},[1 2; 3 4; 1 3; 2 4], 'estimator', 'trimmed mean','newfig','sub');
ylabel('Entropy differences','Fontsize',10); title('Trimmed mean Differences','Fontsize',12); xlabel('')
subplot(3,4,6);
[WMd_diff,WMd_dCI,WMd_p,WMd_alphav,h2] = rst_multicompare(WMd{:,:},[1 2; 3 4; 1 3; 2 4], 'estimator', 'trimmed mean','newfig','sub');
ylabel('Entropy differences','Fontsize',10); title('Trimmed mean Differences','Fontsize',12); xlabel('')
subplot(3,4,10);
[CSFd_diff,CSFd_dCI,CSFd_p,CSFd_alphav,h3] = rst_multicompare(CSFd{:,:},[1 2; 3 4; 1 3; 2 4], 'estimator', 'trimmed mean','newfig','sub');
ylabel('Entropy differences','Fontsize',10); title('Trimmed mean Differences','Fontsize',12); xlabel('')

GMd_table = table(GMd_dCI(:,:)',GMd_p','VariableNames',{'CI','p-value'},'RowName',{'T1nG1 vs T1nG2', 'T12nG1 vs T12nG2', 'T1nG1 vs T12nG1', 'T1nG2 vs T12nG2'});
WMd_table = table(WMd_dCI(:,:)',WMd_p','VariableNames',{'CI','p-value'},'RowName',{'T1nG1 vs T1nG2', 'T12nG1 vs T12nG2', 'T1nG1 vs T12nG1', 'T1nG2 vs T12nG2'});
CSFd_table = table(CSFd_dCI(:,:)',CSFd_p','VariableNames',{'CI','p-value'},'RowName',{'T1nG1 vs T1nG2', 'T12nG1 vs T12nG2', 'T1nG1 vs T12nG1', 'T1nG2 vs T12nG2'});

disp("Entropy")
disp('-------------------------')
disp('     Discovery set')

GMd_change_count  = mean((GMd.T1_nG1-GMd.T1_nG2)<0)*100;
WMd_change_count  = mean((WMd.T1_nG1-WMd.T1_nG2)>0)*100;
CSFd_change_count = mean((CSFd.T1_nG1-CSFd.T1_nG2)<0)*100;
warning("Using T1 only, adding a Gaussian increased GM entropy by %g (%g%% of subjects) and decreased WM entropy by %g (%g%%) and increased CSF entropy by %g (%g%%) (%g increased)", ...
    abs(GMd_est(2)-GMd_est(1)),GMd_change_count,abs(WMd_est(2)-WMd_est(1)),WMd_change_count,abs(CSFd_est(2)-CSFd_est(1)),CSFd_change_count,abs(GMd_est(2)-GMd_est(1)) - abs(WMd_est(2)-WMd_est(1)) + abs(CSFd_est(2)-CSFd_est(1)))

% -------
GMd_change_count  = mean((GMd.T12_nG1-GMd.T12_nG2)>0)*100;
WMd_change_count  = mean((WMd.T12_nG1-WMd.T12_nG2)>0)*100;
CSFd_change_count = mean((CSFd.T12_nG1-CSFd.T12_nG2)>0)*100;
warning('Using T1 and T2, adding a Gaussian decreased GM entropy by %g (%g%% of subjects) and WM entropy by %g (%g%%) and CSF entropy by %g (%g%%) (%g decreased)', ...
    abs(GMd_est(4)-GMd_est(3)),GMd_change_count,abs(WMd_est(4)-WMd_est(3)),WMd_change_count,abs(CSFd_est(4)-CSFd_est(3)),CSFd_change_count,abs(GMd_est(4)-GMd_est(3)) + abs(WMd_est(4)-WMd_est(3)) + abs(CSFd_est(4)-CSFd_est(3)))

% -------
GMd_change_count  = mean((GMd.T1_nG1-GMd.T12_nG1)>0)*100;
WMd_change_count  = mean((WMd.T1_nG1-WMd.T12_nG1)>0)*100;
CSFd_change_count = mean((CSFd.T1_nG1-CSFd.T12_nG1)<0)*100;
warning('With 1 Gaussian only, adding the T2 image decreased GM entropy by %g (%g%% of subjects) and WM entropy by %g (%g%%) and increased CSF entropy by %g (%g%%) (%g increased)', ...
    abs(GMd_est(3)-GMd_est(1)),GMd_change_count,abs(WMd_est(3)-WMd_est(1)),WMd_change_count,abs(CSFd_est(3)-CSFd_est(1)),CSFd_change_count,abs(abs(GMd_est(3)-GMd_est(1)) + abs(WMd_est(3)-WMd_est(1)) - abs(CSFd_est(3)-CSFd_est(1))))

% -------
GMd_change_count  = mean((GMd.T1_nG2-GMd.T12_nG2)>0)*100;
WMd_change_count  = mean((WMd.T1_nG2-WMd.T12_nG2)>0)*100;
CSFd_change_count = mean((CSFd.T1_nG2-CSFd.T12_nG2)>0)*100;
warning('With 2 Gaussians, adding the T2 image decreased GM entropy by %g (%g%% of subjects) and WM entropy by %g (%g%%) and CSF entropy by %g (%g%%) (%g decreased)', ...
    abs(GMd_est(4)-GMd_est(2)),GMd_change_count,abs(WMd_est(4)-WMd_est(2)),WMd_change_count,abs(CSFd_est(4)-CSFd_est(2)),CSFd_change_count,abs(GMd_est(4)-GMd_est(2)) + abs(WMd_est(4)-WMd_est(2)) + abs(CSFd_est(4)-CSFd_est(2)))

disp(DI_table)
disp(GMd_table)
disp(WMd_table)
disp(CSFd_table)
disp("--------")
disp("By adding a second Gaussian")
disp("decreased entropy.")
disp("But by only using one Gaussian an increased the entropy,")
disp("which means that the tissue's voxel probabilities")
disp("are distributed across a wider range of values")
disp("--------")

% replication set
subplot(3,4,3);
[GMt_est, CIt_GM]   = rst_data_plot(GMt{:,:}, 'estimator','trimmed mean','newfig','sub');
title('Grey Matter test set','Fontsize',12); ylabel('GM volumes'); subplot(3,4,7);
[WMt_est, CIt_WM]   = rst_data_plot(WMt{:,:}, 'estimator','trimmed mean','newfig','sub');
title('White Matter test set','Fontsize',12); ylabel('WM volumes'); subplot(3,4,11);
[CSFt_est, CIt_CSF] = rst_data_plot(CSFt{:,:}, 'estimator','trimmed mean','newfig','sub');
title('CSF test set','Fontsize',12); ylabel('CSF volumes');

TrimmedMeans = [GMt_est; WMt_est; CSFt_est];
LowerConfs   = [CIt_GM(1,:); CIt_WM(1,:); CIt_CSF(1,:)];
HigherConfs  = [CIt_GM(2,:); CIt_WM(2,:); CIt_CSF(2,:)];

T1_nG1  = [LowerConfs(:,1) TrimmedMeans(:,1) HigherConfs(:,1)];
T1_nG2  = [LowerConfs(:,2) TrimmedMeans(:,2) HigherConfs(:,2)];
T12_nG1 = [LowerConfs(:,3) TrimmedMeans(:,3) HigherConfs(:,3)];
T12_nG2 = [LowerConfs(:,4) TrimmedMeans(:,4) HigherConfs(:,4)];

DI_table = table(T1_nG1,T1_nG2,T12_nG1,T12_nG2,...
 'RowNames',{'GM','WM','CSF'});

subplot(3,4,4);
Data = [GMt{:,1}-GMt{:,2}, GMt{:,3}-GMt{:,4},...
    GMt{:,1}-GMt{:,3},GMt{:,2}-GMt{:,4}];
[h1,GMt_dCI,GMt_p] = rst_1ttest(Data,'estimator','trimmed mean','newfig','no');
GMt_diff = rst_trimmean(Data);
ylabel('Entropy differences','Fontsize',10); title('Trimmed mean Differences','Fontsize',12);
subplot(3,4,8);
Data = [WMt{:,1}-WMt{:,2}, WMt{:,3}-WMt{:,4},...
    WMt{:,1}-WMt{:,3},WMt{:,2}-WMt{:,4}];
[h2,WMt_dCI,WMt_p] = rst_1ttest(Data,'estimator','trimmed mean','newfig','no');
WMt_diff = rst_trimmean(Data);
ylabel('Entropy differences','Fontsize',10); title('Trimmed mean Differences','Fontsize',12);
subplot(3,4,12);
Data = [CSFt{:,1}-CSFt{:,2}, CSFt{:,3}-CSFt{:,4},...
    CSFt{:,1}-CSFt{:,3},CSFt{:,2}-CSFt{:,4}];
[h3,CSFt_dCI,CSFt_p] = rst_1ttest(Data,'estimator','trimmed mean','newfig','no');
CSFt_diff = rst_trimmean(Data);
ylabel('Entropy differences','Fontsize',10); title('Trimmed mean Differences','Fontsize',12);

GMt_table  = table(GMt_dCI(:,:)',GMt_p','VariableNames',{'CI','p-value'},'RowName',{'T1nG1 vs T1nG2', 'T12nG1 vs T12nG2', 'T1nG1 vs T12nG1', 'T1nG2 vs T12nG2'});
WMt_table  = table(WMt_dCI(:,:)',WMt_p','VariableNames',{'CI','p-value'},'RowName',{'T1nG1 vs T1nG2', 'T12nG1 vs T12nG2', 'T1nG1 vs T12nG1', 'T1nG2 vs T12nG2'});
CSFt_table = table(CSFt_dCI(:,:)',CSFt_p','VariableNames',{'CI','p-value'},'RowName',{'T1nG1 vs T1nG2', 'T12nG1 vs T12nG2', 'T1nG1 vs T12nG1', 'T1nG2 vs T12nG2'});

disp('-------------------------')
disp('    Validation set')

GMt_change_count  = mean((GMt.T1_nG1-GMt.T1_nG2)<0)*100;
WMt_change_count  = mean((WMt.T1_nG1-WMt.T1_nG2)>0)*100;
warning("Using T1 only, adding a Gaussian increased GM entropy by %g (%g%% of subjects) and decreased WM entropy by %g (%g%%), but do not show evidence for changing CSF entropy (%g decreased)", ...
    abs(GMt_est(2)-GMt_est(1)),GMt_change_count,abs(WMt_est(2)-WMt_est(1)),WMt_change_count,abs(abs(GMt_est(2)-GMt_est(1))-abs(WMt_est(2)-WMt_est(1))))

% -------
GMt_change_count  = mean((GMt.T12_nG1-GMt.T12_nG2)>0)*100;
WMt_change_count  = mean((WMt.T12_nG1-WMt.T12_nG2)>0)*100;
CSFt_change_count = mean((CSFt.T12_nG1-CSFt.T12_nG2)>0)*100;
warning('Using T1 and T2, adding a Gaussian decreased GM entropy by %g (%g%% of subjects) and WM entropy by %g (%g%%) and CSF entropy by %g (%g%%) (%g decreased)', ...
    abs(GMt_est(4)-GMt_est(3)),GMt_change_count,abs(WMt_est(4)-WMt_est(3)),WMt_change_count,abs(CSFt_est(4)-CSFt_est(3)),CSFt_change_count,abs(GMt_est(4)-GMt_est(3))+abs(WMt_est(4)-WMt_est(3))+abs(CSFt_est(4)-CSFt_est(3)))

% -------
WMt_change_count  = mean((WMt.T1_nG1-WMt.T12_nG1)<0)*100;
CSFt_change_count = mean((CSFt.T1_nG1-CSFt.T12_nG1)>0)*100;
warning('With 1 Gaussian only, adding the T2 image do not show evidence for changing GM entropy and increased WM entropy by %g (%g%%) and decreased CSF entropy by %g (%g%%) (%g decreased)', ...
    abs(WMt_est(3)-WMt_est(1)),WMt_change_count,abs(CSFt_est(3)-CSFt_est(1)),CSFt_change_count,abs(abs(WMt_est(3)-WMt_est(1))-abs(CSFt_est(3)-CSFt_est(1))))

% -------
GMt_change_count  = mean((GMt.T1_nG2-GMt.T12_nG2)>0)*100;
WMt_change_count  = mean((WMt.T1_nG2-WMt.T12_nG2)>0)*100;
CSFt_change_count = mean((CSFt.T1_nG2-CSFt.T12_nG2)>0)*100;
warning('With 2 Gaussians, adding the T2 image decreased GM entropy by %g (%g%% of subjects) and WM entropy by %g (%g%%) and CSF entropy by %g (%g%%) (%g decreased)', ...
    abs(GMt_est(4)-GMt_est(2)),GMt_change_count,abs(WMt_est(4)-WMt_est(2)),WMt_change_count,abs(CSFt_est(4)-CSFt_est(2)),CSFt_change_count,abs(GMt_est(4)-GMt_est(2))+abs(WMt_est(4)-WMt_est(2))+abs(CSFt_est(4)-CSFt_est(2)))

disp(DI_table)
disp(GMt_table)
disp(WMt_table)
disp(CSFt_table)
disp("--------")
disp("The test set could not replicate the findings of discovery set")
disp("The test set don't show that adding the second Gaussian decreases the entropy")
disp("--------")
disp('-------------------------')

clearvars

% -------------------------------------------------------------------------
%% How much of the missing volumes are now vessels? 
% we investigated where missing volumes are located, thus focusing analyses
% on comparing modality changes only (ie T1w vs, T1w and T2w image inputs in 
% the 1 Gaussian and in the 2 Gaussians models)
% -------------------------------------------------------------------------


% where do we see differences spatially
clear variables
Vessels = spm_read_vols(spm_vol(fullfile(fileparts(pwd),...
    ['code' filesep 'Atlases' filesep 'Vessels' filesep 'rvesselRadius.nii'])));
dataset = {'nrudataset','ds003653'};
names = [1 1; 1 2; 12 1; 12 2];
types = {'mean','var'};

for d=1:2
    for op = 1:4
        file        = fullfile(dataset{d},['mean_modalityT' num2str(names(op,1)) '_NGaussian' num2str(names(op,2)) '.nii.gz']);
        V           = spm_vol(gunzip(file));
        MeanImg{op} = spm_read_vols(V{1});
        spm_unlink(V{1}(1).fname)
        file        = fullfile(dataset{d},['var_modalityT' num2str(names(op,1)) '_NGaussian' num2str(names(op,2)) '.nii.gz']);
        V           = spm_vol(gunzip(file));
        VarImg{op} = spm_read_vols(V{1});
        spm_unlink(V{1}(1).fname)
    end

    new = V{1}(1);
    GMd = readtable([dataset{d} filesep 'GreyMatter_volumes.csv'],'ReadRowNames',false);  % High probability of Grey matter in vessels
    for tissue = 1:3
        T1   = (MeanImg{1}(:,:,:,tissue)+MeanImg{2}(:,:,:,tissue))./2;
        T12  = (MeanImg{3}(:,:,:,tissue)+MeanImg{4}(:,:,:,tissue))./2;
        S21  = (VarImg{1}(:,:,:,tissue)+VarImg{2}(:,:,:,tissue))./2;
        S212 = (VarImg{3}(:,:,:,tissue)+VarImg{4}(:,:,:,tissue))./2;
        D    = T1-T12;
        S    = sqrt((S21+S212)./2); % not ideal, std of data rater than of diff
        T    = D./(S./sqrt(size(GMd,1)));
        P    = 2 * tcdf(-abs(T), size(GMd,1) - 1);
        allp = sort(P(:));
        V    = length(allp);
        I    = (1:V)';
        cVN  = sum(1./(1:V));
        pN   = allp(max(find(allp<=I/V*0.05/cVN))); % FDR threshold
        new.descrip         = ['T-test tissue class ' num2str(tissue)];
        new.private.descrip = '3D';
        new.fname           = [dataset{d} filesep 'T-test_tissue_class_' num2str(tissue) '.nii'];
        spm_write_vol(new,T);
        new.fname           = [dataset{d} filesep 'T-test_tissue_class_' num2str(tissue) 'thresholded.nii'];
        spm_write_vol(new,T.*(P<pN));
        new.fname           = [dataset{d} filesep 'T-test_tissue_class_' num2str(tissue) 'thresholded_masked.nii'];
        spm_write_vol(new,T.*(P<pN).*(Vessels>.5));        
    end
end

%% How much of the missing volumes are now bone (Class 4)?
cd('../')
%num_voxels = NaN(2,4);

% Discovery set
%dartel = spm_read_vols(spm_vol(fullfile('nrudataset', filesep, 'sub-52334', filesep, 'anat', filesep, 'templateT1_nG1_6.nii,4')));  num_voxels(1,1) = nnz(dartel); % Count the non-zero voxels
%dartel = spm_read_vols(spm_vol(fullfile('nrudataset', filesep, 'sub-52334', filesep, 'anat', filesep, 'templateT1_nG2_6.nii,4')));  num_voxels(1,2) = nnz(dartel); % Count the non-zero voxels
%dartel = spm_read_vols(spm_vol(fullfile('nrudataset', filesep, 'sub-52334', filesep, 'anat', filesep, 'templateT12_nG1_6.nii,4'))); num_voxels(1,3) = nnz(dartel); % Count the non-zero voxels
%dartel = spm_read_vols(spm_vol(fullfile('nrudataset', filesep, 'sub-52334', filesep, 'anat', filesep, 'templateT12_nG2_6.nii,4'))); num_voxels(1,4) = nnz(dartel); % Count the non-zero voxels

disp("Comparing dartel template TIV with the calculated TIV for discovery")
disp(spm_summarise(fullfile('nrudataset', filesep, 'sub-52334', filesep, 'anat', filesep, 'templateT1_nG1_6.nii,1'), 'all', 'litres')+ ...
    spm_summarise(fullfile('nrudataset', filesep, 'sub-52334', filesep, 'anat', filesep, 'templateT1_nG1_6.nii,2'), 'all', 'litres')+ ...
    spm_summarise(fullfile('nrudataset', filesep, 'sub-52334', filesep, 'anat', filesep, 'templateT1_nG1_6.nii,3'), 'all', 'litres'))

disp("Comparing dartel template TIV with the calculated TIV for validation")
disp(spm_summarise(fullfile('ds003653', filesep, 'sub-718216', filesep, 'anat', filesep, 'templateT1_nG1_6.nii,1'), 'all', 'litres')+ ...
    spm_summarise(fullfile('ds003653', filesep, 'sub-718216', filesep, 'anat', filesep, 'templateT1_nG1_6.nii,2'), 'all', 'litres')+ ...
    spm_summarise(fullfile('ds003653', filesep, 'sub-718216', filesep, 'anat', filesep, 'templateT1_nG1_6.nii,3'), 'all', 'litres'))

T = table(spm_summarise(fullfile('nrudataset', filesep, 'sub-52334', filesep, 'anat', filesep, 'templateT1_nG1_6.nii,4'), 'all', 'litres')*1000, ...
          spm_summarise(fullfile('nrudataset', filesep, 'sub-52334', filesep, 'anat', filesep, 'templateT1_nG2_6.nii,4'), 'all', 'litres')*1000, ...
          spm_summarise(fullfile('nrudataset', filesep, 'sub-52334', filesep, 'anat', filesep, 'templateT12_nG1_6.nii,4'), 'all', 'litres')*1000, ...
          spm_summarise(fullfile('nrudataset', filesep, 'sub-52334', filesep, 'anat', filesep, 'templateT12_nG2_6.nii,4'), 'all', 'litres')*1000, ...
          'RowName',{'Volume (ml)'}, 'VariableNames', {'T1_nG1', 'T1_nG2', 'T12_nG1', 'T12_nG2'});
      
disp('Discovery set')
disp('-----------------------')
disp(T)
warning('By adding T2w image the tissue class 4 grew with an avage of %g ml',(T.T12_nG1-T.T1_nG1+T.T12_nG2-T.T1_nG2)/2)

% Validation set
%dartel = spm_read_vols(spm_vol(fullfile('ds003653', filesep, 'sub-718216', filesep, 'anat', filesep, 'templateT1_nG1_6.nii,4'))); num_voxels(2,1) = nnz(dartel); % Count the non-zero voxels
%dartel = spm_read_vols(spm_vol(fullfile('ds003653', filesep, 'sub-718216', filesep, 'anat', filesep, 'templateT1_nG2_6.nii,4'))); num_voxels(2,2) = nnz(dartel); % Count the non-zero voxels
%dartel = spm_read_vols(spm_vol(fullfile('ds003653', filesep, 'sub-718216', filesep, 'anat', filesep, 'templateT12_nG1_6.nii,4')));num_voxels(2,3) = nnz(dartel); % Count the non-zero voxels
%dartel = spm_read_vols(spm_vol(fullfile('ds003653', filesep, 'sub-718216', filesep, 'anat', filesep, 'templateT12_nG2_6.nii,4')));num_voxels(2,4) = nnz(dartel); % Count the non-zero voxels

T = table(spm_summarise(fullfile('ds003653', filesep, 'sub-718216', filesep, 'anat', filesep, 'templateT1_nG1_6.nii,4'), 'all', 'litres')*1000, ...
          spm_summarise(fullfile('ds003653', filesep, 'sub-718216', filesep, 'anat', filesep, 'templateT1_nG2_6.nii,4'), 'all', 'litres')*1000, ...
          spm_summarise(fullfile('ds003653', filesep, 'sub-718216', filesep, 'anat', filesep, 'templateT12_nG1_6.nii,4'), 'all', 'litres')*1000, ...
          spm_summarise(fullfile('ds003653', filesep, 'sub-718216', filesep, 'anat', filesep, 'templateT12_nG2_6.nii,4'), 'all', 'litres')*1000, ...
          'RowName',{'Volume (L)'}, 'VariableNames', {'T1_nG1', 'T1_nG2', 'T12_nG1', 'T12_nG2'});
      
disp('Validation set')
disp('-----------------------')
disp(T)
warning('By adding T2w image the tissue class 4 grew with an avage of %g ml',(T.T12_nG1-T.T1_nG1+T.T12_nG2-T.T1_nG2)/2)
