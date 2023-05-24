% statistical analysis of multispectral segmentation 
% variables are appended with d for discovery (NRU data)
% or with a t for test or validation set (ds003653).

cd('../results')

%% what is the total intracranial volume (TIV) for the four types of segmentation
% -------------------------------------------------------------------------------
GMd  = readtable(['NRU_dataset' filesep 'GrayMatter_volumes.csv'],'ReadRowNames',false);           
WMd  = readtable(['NRU_dataset' filesep 'WhiteMatter_volumes.csv'],'ReadRowNames',false);           
CSFd = readtable(['NRU_dataset' filesep 'CSF_volumes.csv'],'ReadRowNames',false);           
TIVd = [GMd{:,1}+WMd{:,1}+CSFd{:,1} GMd{:,2}+WMd{:,2}+CSFd{:,2} ...
    GMd{:,3}+WMd{:,3}+CSFd{:,3} GMd{:,4}+WMd{:,4}+CSFd{:,4}];
GMt  = readtable(['ds003653' filesep 'GrayMatter_volumes.csv'],'ReadRowNames',false);           
WMt  = readtable(['ds003653' filesep 'WhiteMatter_volumes.csv'],'ReadRowNames',false);           
CSFt = readtable(['ds003653' filesep 'CSF_volumes.csv'],'ReadRowNames',false);           
TIVt = [GMt{:,1}+WMt{:,1}+CSFt{:,1} GMt{:,2}+WMt{:,2}+CSFt{:,2} ...
    GMt{:,3}+WMt{:,3}+CSFt{:,3} GMt{:,4}+WMt{:,4}+CSFt{:,4}];

% in the discovery set test main effects and interaction using a Hotelling
% test (repeated measure ANOVA) and multiple pair differences (alphav is adjusted
% using Hochberg step-up procedure)
result = rst_rep_anova_T2(TIVd,[],[2 2],1000,{'modality','n_gaussians'});
disp('-----');
warning('significant effect of modality, nb of gaussians but no interaction') 
disp(result)
disp('-----')

figure('Name','TIV'); subplot(2,2,1);
[TIVd_est, TIVd_CI]   = rst_data_plot(TIVd, 'estimator','trimmed mean','newfig','sub');
title('TIV discovery set','Fontsize',12); ylabel('volumes'); subplot(2,2,2)
[TIVd_diff,TIVd_dCI,TIVd_p,TIVd_alphav,h] = rst_multicompare(TIVd,[1 2; 3 4; 1 3; 2 4], 'estimator', 'trimmed mean','newfig','sub');
ylabel('volume differences','Fontsize',10); title('Trimmed mean Differences','Fontsize',12); xlabel('')
disp('-----');
warning('adding 1 Gaussian increases TIV by %g and %g ml for unimodal and multimodal segmentation',abs(TIVd_diff(1)),abs(TIVd_diff(2)));
warning('adding a T2 image decreases TIV by %g and %g ml for 1 Gaussian and 2 Gaussians models',TIVd_diff(3),TIVd_diff(4));

% replication set - test for the same differences found as above using
% strict Bonferonni correction
subplot(2,2,3);
[TIVt_est, TIVt_CI]   = rst_data_plot(TIVt, 'estimator','trimmed mean','newfig','sub');
title('TIV validation set','Fontsize',12); ylabel('volumes'); 
xlabel('T1-1G T1-2G T12-1G T12-2G'); subplot(2,2,4)
Data = [TIVt(:,1)-TIVt(:,2), TIVt(:,3)-TIVt(:,4),...
    TIVt(:,1)-TIVt(:,3),TIVt(:,2)-TIVt(:,4)];
[h,CI,p] = rst_1ttest(Data,'estimator','trimmed mean','newfig','no');
TIVt_diff = rst_trimmean(Data);
ylabel('volume differences','Fontsize',10); title('Trimmed mean Differences','Fontsize',12);
xlabel('1vs.2 3vs.4 1vs.3 2vs.4')
disp('-----');
warning('validation set confirms differences observed in the discovery set');
warning('adding 1 Gaussian increases TIV by %g and %g ml for unimodal and multimodal segmentation',abs(TIVt_diff(1)),abs(TIVt_diff(2)));
warning('adding a T2 image decreases TIV by %g and %g ml for 1 Gaussian and 2 Gaussians models',TIVt_diff(3),TIVt_diff(4));

% summary table
summary = table([TIVd_CI(1,1) TIVd_est(1) TIVd_CI(2,1); TIVt_CI(1,1) TIVt_est(1) TIVt_CI(2,1)],...
    [TIVd_CI(1,2) TIVd_est(2) TIVd_CI(2,2); TIVt_CI(1,2) TIVt_est(2) TIVt_CI(2,2)],...
    [TIVd_CI(1,3) TIVd_est(3) TIVd_CI(2,3); TIVt_CI(1,3) TIVt_est(3) TIVt_CI(2,3)],...
    [TIVd_CI(1,4) TIVd_est(4) TIVd_CI(2,4); TIVt_CI(1,4) TIVt_est(4) TIVt_CI(2,4)],...
    'RowNames',{'discovery','test'},'VariableNames',{'T1_nG1','T1_nG2','T12_nG1','T12_nG2'});
disp(summary); warning('note how with 1 Gaussian distributions barely overlap')

%% how are changes in TIV explained by tissue type
% ------------------------------------------------

figure('Name','Tissue volumes'); 
subplot(3,4,1);
[GMd_est, CId_GM]   = rst_data_plot(GMd{:,:}, 'estimator','trimmed mean','newfig','sub');
title('Gray Matter discovery set','Fontsize',12); ylabel('GM volumes'); subplot(3,4,5);
[WMd_est, CId_WM]   = rst_data_plot(WMd{:,:}, 'estimator','trimmed mean','newfig','sub');
title('White Matter discovery set','Fontsize',12); ylabel('WM volumes'); subplot(3,4,9);
[CSFd_est, CId_CSF] = rst_data_plot(CSFd{:,:}, 'estimator','trimmed mean','newfig','sub');
title('CSF discovery set','Fontsize',12); ylabel('CSF volumes');

subplot(3,4,2);
[GMd_diff,GMd_dCI,GMd_p,GMd_alphav,h1] = rst_multicompare(GMd{:,:},[1 2; 3 4; 1 3; 2 4], 'estimator', 'trimmed mean','newfig','sub');
ylabel('volume differences','Fontsize',10); title('Trimmed mean Differences','Fontsize',12); xlabel('')
subplot(3,4,6);
[WMd_diff,WMd_dCI,WMd_p,WMd_alphav,h2] = rst_multicompare(WMd{:,:},[1 2; 3 4; 1 3; 2 4], 'estimator', 'trimmed mean','newfig','sub');
ylabel('volume differences','Fontsize',10); title('Trimmed mean Differences','Fontsize',12); xlabel('')
subplot(3,4,10);
[CSFd_diff,CSFd_dCI,CSFd_p,CSFd_alphav,h3] = rst_multicompare(CSFd{:,:},[1 2; 3 4; 1 3; 2 4], 'estimator', 'trimmed mean','newfig','sub');
ylabel('volume differences','Fontsize',10); title('Trimmed mean Differences','Fontsize',12); xlabel('')

% look at the relationship between tissue classes using correlations, what
% does it tell us about how changing model input or parameters changes
% tissue class attributions
[rd,td,pvald,rCId,alphav] = Spearman([GMd{:,:} GMd{:,:} WMd{:,:}],[WMd{:,:} CSFd{:,:} CSFd{:,:}]);
[rt,tt,pvalt,rCIt,alphav] = Spearman([GMt{:,:} GMt{:,:} WMt{:,:}],[WMt{:,:} CSFt{:,:} CSFt{:,:}]);

% results
% -------
GMd_change_count  = mean((GMd.T1_nG1-GMd.T1_nG2)>0)*100;
WMd_change_count  = mean((WMd.T1_nG1-WMd.T1_nG2)<0)*100;
CSFd_change_count = mean((CSFd.T1_nG1-CSFd.T1_nG2)<0)*100;
warning('Using T1 only, adding a Gaussian decreased GM volumes by %g ml (%g%% of subjects) and increase WM by %g ml (%g%%) and CSF by %g ml (%g%%) changes in proportion (diff = %g ml)', ...
    GMd_diff(1),GMd_change_count,abs(WMd_diff(1)),WMd_change_count,abs(CSFd_diff(1)),CSFd_change_count,GMd_diff(1) - abs(WMd_diff(1)) - abs(CSFd_diff(1)))
warning('the relationships between tissue classes are unchanged')

% -------
GMd_change_count  = mean((GMd.T12_nG1-GMd.T12_nG2)<0)*100;
WMd_change_count  = mean((WMd.T12_nG1-WMd.T12_nG2)<0)*100;
CSFd_change_count = mean((CSFd.T12_nG1-CSFd.T12_nG2)>0)*100;
warning('Using T1 and T2, adding a Gaussian increased GM by %g ml (%g%% of subjects) and WM by %g ml (%g%%) but decreased CSF by %g ml (%g%%) in proportion (diff = %g ml)', ...
    abs(GMd_diff(2)),GMd_change_count,abs(WMd_diff(2)),WMd_change_count,abs(CSFd_diff(2)),CSFd_change_count, abs(GMd_diff(2)) + abs(WMd_diff(2)) - abs(CSFd_diff(2)))

% -------
GMd_change_count  = mean((GMd.T1_nG1-GMd.T12_nG1)>0)*100;
WMd_change_count  = mean((WMd.T1_nG1-WMd.T12_nG1)>0)*100;
CSFd_change_count = mean((CSFd.T1_nG1-CSFd.T12_nG1)<0)*100;
warning('With 1 Gaussian only, adding the T2 image decreased GM by %g ml (%g%% of subjects) and WM by %g ml (%g%%) but increases CSF by %g ml (%g%%) and unproportionally (%g ml missing)', ...
    GMd_diff(3),GMd_change_count,WMd_diff(3),WMd_change_count,abs(CSFd_diff(3)),CSFd_change_count,GMd_diff(3) + WMd_diff(3) - abs(CSFd_diff(3)))

% -------
GMd_change_count  = mean((GMd.T1_nG2-GMd.T12_nG2)<0)*100;
WMd_change_count  = mean((WMd.T1_nG2-WMd.T12_nG2)>0)*100;
CSFd_change_count = mean((CSFd.T1_nG2-CSFd.T12_nG2)>0)*100;
warning('With 2 Gaussians, adding the T2 image increased GM volume by %g ml (%g%% of subjects) and decreased WM by %g ml (%g%%) and CSF by %g ml (%g%%) and unproportionally (%g ml missing)', ...
    abs(GMd_diff(4)),GMd_change_count,WMd_diff(4),WMd_change_count,abs(CSFd_diff(4)),CSFd_change_count, WMd_diff(4) + CSFd_diff(4) - abs(GMd_diff(4)))

figure('Name','Adding a Gaussian to T1 input'); set(gcf,'Color','w'); 

subplot(4,3,1);
scatter(GMd.T1_nG1,WMd.T1_nG1,100,'filled'); grid on; box on; xlabel('X','FontSize',10); ylabel('Y','FontSize',10);
M = sprintf('GM/WM T1 nG1 r=%g \n %g%%CI [%.2f %.2f]',rd(1),(1-alphav)*100,rCId(1,1),rCId(2,1));
title(M,'FontSize',12); h=lsline; set(h,'Color','r','LineWidth',4);
subplot(4,3,4);
scatter(GMd.T1_nG2,WMd.T1_nG2,100,'filled'); grid on; box on; xlabel('X','FontSize',10); ylabel('Y','FontSize',10);
M = sprintf('GM/WM T1 nG2 r=%g \n %g%%CI [%.2f %.2f]',rd(2),(1-alphav)*100,rCId(1,2),rCId(2,2));
title(M,'FontSize',12); h=lsline; set(h,'Color','r','LineWidth',4);
subplot(4,3,7);
scatter(GMd.T12_nG1,WMd.T12_nG1,100,'filled'); grid on; box on; xlabel('X','FontSize',10); ylabel('Y','FontSize',10);
M = sprintf('GM/WM T1-T2 nG1 r=%g \n %g%%CI [%.2f %.2f]',rd(3),(1-alphav)*100,rCId(1,3),rCId(2,3));
title(M,'FontSize',12); h=lsline; set(h,'Color','r','LineWidth',4);
subplot(4,3,10);
scatter(GMd.T12_nG2,WMd.T12_nG2,100,'filled'); grid on; box on; xlabel('X','FontSize',10); ylabel('Y','FontSize',10);
M = sprintf('GM/WM T1-T2 nG2 r=%g \n %g%%CI [%.2f %.2f]',rd(4),(1-alphav)*100,rCId(1,4),rCId(2,4));
title(M,'FontSize',12); h=lsline; set(h,'Color','r','LineWidth',4);

subplot(4,3,2);
scatter(GMd.T1_nG1,CSFd.T1_nG1,100,'filled'); grid on; box on; xlabel('X','FontSize',10); ylabel('Y','FontSize',10);
M = sprintf('GM/CSF T1 nG1 r=%g \n %g%%CI [%.2f %.2f]',rd(5),(1-alphav)*100,rCId(1,5),rCId(2,5));
title(M,'FontSize',12); h=lsline; set(h,'Color','r','LineWidth',4);
subplot(4,3,5);
scatter(GMd.T1_nG2,CSFd.T1_nG2,100,'filled'); grid on; box on; xlabel('X','FontSize',10); ylabel('Y','FontSize',10);
M = sprintf('GM/CSF T1 nG2 r=%g \n %g%%CI [%.2f %.2f]',rd(6),(1-alphav)*100,rCId(1,6),rCId(2,6));
title(M,'FontSize',12); h=lsline; set(h,'Color','r','LineWidth',4);
subplot(4,3,8);
scatter(GMd.T12_nG1,CSFd.T12_nG1,100,'filled'); grid on; box on; xlabel('X','FontSize',10); ylabel('Y','FontSize',10);
M = sprintf('GM/CSF T1-T2 nG1 r=%g \n %g%%CI [%.2f %.2f]',rd(7),(1-alphav)*100,rCId(1,7),rCId(2,7));
title(M,'FontSize',12); h=lsline; set(h,'Color','r','LineWidth',4);
subplot(4,3,11);
scatter(GMd.T12_nG2,CSFd.T12_nG2,100,'filled'); grid on; box on; xlabel('X','FontSize',10); ylabel('Y','FontSize',10);
M = sprintf('GM/CSF T1-T2 nG2 r=%g \n %g%%CI [%.2f %.2f]',rd(8),(1-alphav)*100,rCId(1,8),rCId(2,8));
title(M,'FontSize',12); h=lsline; set(h,'Color','r','LineWidth',4);

subplot(4,3,3);
scatter(WMd.T1_nG1,CSFd.T1_nG1,100,'filled'); grid on; box on; xlabel('X','FontSize',10); ylabel('Y','FontSize',10);
M = sprintf('WM/CSF T1 nG1 r=%g \n %g%%CI [%.2f %.2f]',rd(9),(1-alphav)*100,rCId(1,9),rCId(2,9));
title(M,'FontSize',12); h=lsline; set(h,'Color','r','LineWidth',4);
subplot(4,3,6);
scatter(WMd.T1_nG2,CSFd.T1_nG2,100,'filled'); grid on; box on; xlabel('X','FontSize',10); ylabel('Y','FontSize',10);
M = sprintf('WM/CSF T1 nG2 r=%g \n %g%%CI [%.2f %.2f]',rd(10),(1-alphav)*100,rCId(1,10),rCId(2,10));
title(M,'FontSize',12); h=lsline; set(h,'Color','r','LineWidth',4);
subplot(4,3,9);
scatter(WMd.T12_nG1,CSFd.T12_nG1,100,'filled'); grid on; box on; xlabel('X','FontSize',10); ylabel('Y','FontSize',10);
M = sprintf('WM/CSF T1-T2 nG1 r=%g \n %g%%CI [%.2f %.2f]',rd(11),(1-alphav)*100,rCId(1,11),rCId(2,11));
title(M,'FontSize',12); h=lsline; set(h,'Color','r','LineWidth',4);
subplot(4,3,12);
scatter(WMd.T12_nG2,CSFd.T12_nG2,100,'filled'); grid on; box on; xlabel('X','FontSize',10); ylabel('Y','FontSize',10);
M = sprintf('WM/CSF T1-T2 nG2 r=%g \n %g%%CI [%.2f %.2f]',rd(12),(1-alphav)*100,rCId(1,12),rCId(2,12));
title(M,'FontSize',12); h=lsline; set(h,'Color','r','LineWidth',4);

disp('------------------')
disp('the proportional changes in tissue volumes adding a Gaussian while')
disp('the correlation structure is preserved, indicates a change in voxel')
disp('distribution fitting (mixture) - the unproportional changes in tisue volumes')
disp('and dispruption of correlations indicates a change in tissue class belonging')
disp('with some tissue missing!')
disp('------------------')

% validation set
% -----------------
figure(findobj( 'Type', 'Figure', 'Name', 'Tissue volumes' ));
subplot(3,4,3);
[GMt_est, CIt_GM]   = rst_data_plot(GMt{:,:}, 'estimator','trimmed mean','newfig','sub');
title('Gray Matter test set','Fontsize',12); ylabel('GM volumes'); subplot(3,4,7);
[WMt_est, CIt_WM]   = rst_data_plot(WMt{:,:}, 'estimator','trimmed mean','newfig','sub');
title('White Matter test set','Fontsize',12); ylabel('WM volumes'); subplot(3,4,11);
[CSFt_est, CIt_CSF] = rst_data_plot(CSFt{:,:}, 'estimator','trimmed mean','newfig','sub');
title('CSF test set','Fontsize',12); ylabel('CSF volumes');

subplot(3,4,4);
Data = [GMt{:,1}-GMt{:,2}, GMt{:,3}-GMt{:,4},...
    GMt{:,1}-GMt{:,3},GMt{:,2}-GMt{:,4}];
[h1,CI,p] = rst_1ttest(Data,'estimator','trimmed mean','newfig','no');
GMt_diff = rst_trimmean(Data);
ylabel('volume differences','Fontsize',10); title('Trimmed mean Differences','Fontsize',12);
subplot(3,4,8);
Data = [WMt{:,1}-WMt{:,2}, WMt{:,3}-WMt{:,4},...
    WMt{:,1}-WMt{:,3},WMt{:,2}-WMt{:,4}];
[h2,CI,p] = rst_1ttest(Data,'estimator','trimmed mean','newfig','no');
WMt_diff = rst_trimmean(Data);
ylabel('volume differences','Fontsize',10); title('Trimmed mean Differences','Fontsize',12);
subplot(3,4,12);
Data = [CSFt{:,1}-CSFt{:,2}, CSFt{:,3}-CSFt{:,4},...
    CSFt{:,1}-CSFt{:,3},CSFt{:,2}-CSFt{:,4}];
[h3,CI,p] = rst_1ttest(Data,'estimator','trimmed mean','newfig','no');
CSFt_diff = rst_trimmean(Data);
ylabel('volume differences','Fontsize',10); title('Trimmed mean Differences','Fontsize',12);

% --------------
GMt_change_count  = mean((GMt.T1_nG1-GMt.T1_nG2)>0)*100;
WMt_change_count  = mean((WMt.T1_nG1-WMt.T1_nG2)<0)*100;
CSFt_change_count = mean((CSFt.T1_nG1-CSFt.T1_nG2)<0)*100;
warning('Using T1 only, adding a Gaussian decreased GM volumes by %g ml (%g%% of subjects) and increase WM by %g ml (%g%%) and CSF by %g ml (%g%%) changes in proportion (diff = %g ml)', ...
    GMt_diff(1),GMt_change_count,abs(WMt_diff(1)),WMt_change_count,abs(CSFt_diff(1)),CSFt_change_count,GMt_diff(1) - abs(WMt_diff(1)) - abs(CSFt_diff(1)))

% -------
GMt_change_count  = mean((GMt.T12_nG1-GMt.T12_nG2)<0)*100;
WMt_change_count  = mean((WMt.T12_nG1-WMt.T12_nG2)<0)*100;
CSFt_change_count = mean((CSFt.T12_nG1-CSFt.T12_nG2)>0)*100;
warning('Using T1 and T2, adding a Gaussian increased GM by %g ml (%g%% of subjects) and WM by %g ml (%g%%) but decreased CSF by %g ml (%g%%) in proportion (diff = %g ml)', ...
    abs(GMt_diff(2)),GMt_change_count,abs(WMt_diff(2)),WMt_change_count,abs(CSFt_diff(2)),CSFt_change_count, abs(GMt_diff(2)) + abs(WMt_diff(2)) - abs(CSFt_diff(2)))

% -------
GMt_change_count  = mean((GMt.T1_nG1-GMt.T12_nG1)>0)*100;
WMt_change_count  = mean((WMt.T1_nG1-WMt.T12_nG1)>0)*100;
CSFt_change_count = mean((CSFt.T1_nG1-CSFt.T12_nG1)>0)*100;
warning('With 1 Gaussian only, adding the T2 image decreased GM by %g ml (%g%% of subjects) and WM by %g ml (%g%%) but also CSF by %g ml (%g%%) and unproportionally (%g ml missing)', ...
    GMt_diff(3),GMt_change_count,WMt_diff(3),WMt_change_count,CSFt_diff(3),CSFt_change_count,GMt_diff(3) + WMt_diff(3) + CSFt_diff(3))

% -------
GMt_change_count  = mean((GMt.T1_nG2-GMt.T12_nG2)<0)*100;
WMt_change_count  = mean((WMt.T1_nG2-WMt.T12_nG2)>0)*100;
CSFt_change_count = mean((CSFt.T1_nG2-CSFt.T12_nG2)>0)*100;
warning('With 2 Gaussians, adding the T2 image hardly changed GM (+%g ml - %g%% of subjects) and WM volumes (-%g ml - %g%%) but decreased CSF by %g ml (%g%%) and unproportionally (%g ml missing)', ...
    abs(GMt_diff(4)),GMt_change_count,abs(WMt_diff(4)),WMt_change_count,abs(CSFt_diff(4)),CSFt_change_count, WMt_diff(4) + CSFt_diff(4) + GMt_diff(4))

% -------------------------------------------------------------------------
%% What does this all mean in terms of similarlity/differences among tissues
% -------------------------------------------------------------------------

% Dunn index
GMd  = readtable(['NRU_dataset' filesep 'GrayMatter_DunnIndexes.csv'],'ReadRowNames',false);           
WMd  = readtable(['NRU_dataset' filesep 'WhiteMatter_DunnIndexes.csv'],'ReadRowNames',false);           
CSFd = readtable(['NRU_dataset' filesep 'CSF_DunnIndexes.csv'],'ReadRowNames',false);           

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

DI_TM_table = table(T1_nG1,T1_nG2,T12_nG1,T12_nG2,...
 'RowNames',{'GM','WM','CSF'});

figure('Name','Dunn Index Trimmed mean'); subplot(3,2,1);
[GMd_diff,GMd_dCI,GMd_p,GMd_alphav,h1] = rst_multicompare(GMd{:,:},[1 2; 3 4; 1 3; 2 4], 'estimator', 'trimmed mean','newfig','sub');
ylabel('Dunn Index differences','Fontsize',10); title('Trimmed mean Differences','Fontsize',12); xlabel('');  subplot(3,2,3);
[WMd_diff,WMd_dCI,WMd_p,WMd_alphav,h2] = rst_multicompare(WMd{:,:},[1 2; 3 4; 1 3; 2 4], 'estimator', 'trimmed mean','newfig','sub');
ylabel('Dunn Index differences','Fontsize',10); title('Trimmed mean Differences','Fontsize',12); xlabel('');  subplot(3,2,5);
[CSFd_diff,CSFd_dCI,CSFd_p,CSFd_alphav,h3] = rst_multicompare(CSFd{:,:},[1 2; 3 4; 1 3; 2 4], 'estimator', 'trimmed mean','newfig','sub');
ylabel('Dunn Index differences','Fontsize',10); title('Trimmed mean Differences','Fontsize',12); xlabel('')

GMd_table = table(GMd_dCI(:,:)',GMd_p','VariableNames',{'CI','p-value'},'RowName',{'T1nG1 vs T1nG2', 'T12nG1 vs T12nG2', 'T1nG1 vs T12nG1', 'T1nG2 vs T12nG2'});
WMd_table = table(WMd_dCI(:,:)',WMd_p','VariableNames',{'CI','p-value'},'RowName',{'T1nG1 vs T1nG2', 'T12nG1 vs T12nG2', 'T1nG1 vs T12nG1', 'T1nG2 vs T12nG2'});
CSFd_table = table(CSFd_dCI(:,:)',CSFd_p','VariableNames',{'CI','p-value'},'RowName',{'T1nG1 vs T1nG2', 'T12nG1 vs T12nG2', 'T1nG1 vs T12nG1', 'T1nG2 vs T12nG2'});

disp("Dunn Index")
disp('-------------------------')
disp('     Discovery set')

GMd_change_count  = mean((GMd.T1_nG1-GMd.T1_nG2)==0)*100;
WMd_change_count  = mean((WMd.T1_nG1-WMd.T1_nG2)==0)*100;
CSFd_change_count = mean((CSFd.T1_nG1-CSFd.T1_nG2)==0)*100;
warning("Using T1 only, adding a Gaussian do not change GM Dunn Index (%g%% of subjects) and WM Dunn Index (%g%%) and CSF Dunn Index (%g%%)", ...
    GMd_change_count,WMd_change_count,CSFd_change_count)

% -------
GMd_change_count  = mean((GMd.T12_nG1-GMd.T12_nG2)==0)*100;
WMd_change_count  = mean((WMd.T12_nG1-WMd.T12_nG2)==0)*100;
CSFd_change_count = mean((CSFd.T12_nG1-CSFd.T12_nG2)==0)*100;
warning('Using T1 and T2, adding a Gaussian do not change GM Dunn Index (%g%% of subjects) and WM Dunn Index (%g%%) and CSF Dunn Index (%g%%)', ...
    GMd_change_count,WMd_change_count,CSFd_change_count)

% -------
GMd_change_count  = mean((GMd.T1_nG1-GMd.T12_nG1)==0)*100;
WMd_change_count  = mean((WMd.T1_nG1-WMd.T12_nG1)==0)*100;
CSFd_change_count = mean((CSFd.T1_nG1-CSFd.T12_nG1)==0)*100;
warning('With 1 Gaussian only, adding the T2 image do not change GM Dunn Index (%g%% of subjects) and WM Dunn Index (%g%%) and CSF Dunn Index (%g%%)', ...
    GMd_change_count,WMd_change_count,CSFd_change_count)

% -------
GMd_change_count  = mean((GMd.T1_nG2-GMd.T12_nG2)==0)*100;
WMd_change_count  = mean((WMd.T1_nG2-WMd.T12_nG2)==0)*100;
CSFd_change_count = mean((CSFd.T1_nG2-CSFd.T12_nG2)==0)*100;
warning('With 2 Gaussians, adding the T2 image do not change GM Dunn Index (%g%% of subjects) and WM Dunn Index (%g%%) and CSF Dunn Index (%g%%)', ...
    GMd_change_count,WMd_change_count,CSFd_change_count)

disp(DI_TM_table)
disp(GMd_table)
disp(WMd_table)
disp(CSFd_table)

disp("--------")
disp("There are not any evidence for")
disp("the conditions changing the Dunn Index")
disp("This tell that the tissues dont drastic changes the ranges")
disp("between the different conditions")
disp("--------")

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

subplot(3,2,2);
Data = [GMt{:,1}-GMt{:,2}, GMt{:,3}-GMt{:,4},...
    GMt{:,1}-GMt{:,3},GMt{:,2}-GMt{:,4}];
[h1,GMt_CI,GMt_p] = rst_1ttest(Data,'estimator','trimmed mean','newfig','no');
GMt_diff = rst_trimmean(Data);
ylabel('Dunn Index differences','Fontsize',10); title('Trimmed mean Differences','Fontsize',12);
subplot(3,2,4);
Data = [WMt{:,1}-WMt{:,2}, WMt{:,3}-WMt{:,4},...
    WMt{:,1}-WMt{:,3},WMt{:,2}-WMt{:,4}];
[h2,WMt_CI,WMt_p] = rst_1ttest(Data,'estimator','trimmed mean','newfig','no');
WMt_diff = rst_trimmean(Data);
ylabel('Dunn Index differences','Fontsize',10); title('Trimmed mean Differences','Fontsize',12);
subplot(3,2,6);
Data = [CSFt{:,1}-CSFt{:,2}, CSFt{:,3}-CSFt{:,4},...
    CSFt{:,1}-CSFt{:,3},CSFt{:,2}-CSFt{:,4}];
[h3,CSFt_CI,CSFt_p] = rst_1ttest(Data,'estimator','trimmed mean','newfig','no');
CSFt_diff = rst_trimmean(Data);
ylabel('Dunn Index differences','Fontsize',10); title('Trimmed mean Differences','Fontsize',12);

GMt_table = table(GMt_CI(:,:)',GMt_p','VariableNames',{'CI','p-value'},'RowName',{'T1nG1 vs T1nG2', 'T12nG1 vs T12nG2', 'T1nG1 vs T12nG1', 'T1nG2 vs T12nG2'});
WMt_table = table(WMt_CI(:,:)',WMt_p','VariableNames',{'CI','p-value'},'RowName',{'T1nG1 vs T1nG2', 'T12nG1 vs T12nG2', 'T1nG1 vs T12nG1', 'T1nG2 vs T12nG2'});
CSFt_table = table(CSFt_CI(:,:)',CSFt_p','VariableNames',{'CI','p-value'},'RowName',{'T1nG1 vs T1nG2', 'T12nG1 vs T12nG2', 'T1nG1 vs T12nG1', 'T1nG2 vs T12nG2'});

disp('-------------------------')
disp('     Validation set')

GMt_change_count  = mean((GMt.T1_nG1-GMt.T1_nG2)==0)*100;
WMt_change_count  = mean((WMt.T1_nG1-WMt.T1_nG2)==0)*100;
CSFt_change_count = mean((CSFt.T1_nG1-CSFt.T1_nG2)==0)*100;
warning("Using T1 only, adding a Gaussian do not change GM Dunn Index (%g%% of subjects) and WM Dunn Index (%g%%) and CSF Dunn Index (%g%%)", ...
    GMt_change_count,WMt_change_count,CSFt_change_count)

% -------
GMt_change_count  = mean((GMt.T12_nG1-GMt.T12_nG2)==0)*100;
WMt_change_count  = mean((WMt.T12_nG1-WMt.T12_nG2)==0)*100;
CSFt_change_count = mean((CSFt.T12_nG1-CSFt.T12_nG2)==0)*100;
warning('Using T1 and T2, adding a Gaussian do not change GM Dunn Index (%g%% of subjects) and WM Dunn Index (%g%%) and CSF Dunn Index (%g%%)', ...
    GMt_change_count,WMt_change_count,CSFt_change_count)

% -------
GMt_change_count  = mean((GMt.T1_nG1-GMt.T12_nG1)==0)*100;
WMt_change_count  = mean((WMt.T1_nG1-WMt.T12_nG1)==0)*100;
CSFt_change_count = mean((CSFt.T1_nG1-CSFt.T12_nG1)==0)*100;
warning('With 1 Gaussian only, adding the T2 image do not change GM Dunn Index (%g%% of subjects) and WM Dunn Index (%g%%) and CSF Dunn Index (%g%%)', ...
    GMt_change_count,WMt_change_count,CSFt_change_count)

% -------
GMt_change_count  = mean((GMt.T1_nG2-GMt.T12_nG2)==0)*100;
WMt_change_count  = mean((WMt.T1_nG2-WMt.T12_nG2)==0)*100;
CSFt_change_count = mean((CSFt.T1_nG2-CSFt.T12_nG2)==0)*100;
warning('With 2 Gaussians, adding the T2 image do not change GM Dunn Index (%g%% of subjects) and WM Dunn Index (%g%%) and CSF Dunn Index (%g%%)', ...
    GMt_change_count,WMt_change_count,CSFt_change_count)

disp(DI_table)
disp(GMt_table)
disp(WMt_table)
disp(CSFt_table)
disp("--------")
disp("Validation set confirmes the findings in the discovery set")
disp("--------")
disp('-------------------------')

clearvars

% Entropy
GMd  = readtable(['NRU_dataset' filesep 'GrayMatter_entropy.csv'],'ReadRowNames',false);           
WMd  = readtable(['NRU_dataset' filesep 'WhiteMatter_entropy.csv'],'ReadRowNames',false);           
CSFd = readtable(['NRU_dataset' filesep 'CSF_entropy.csv'],'ReadRowNames',false);           

GMt  = readtable(['ds003653' filesep 'GrayMatter_entropy.csv'],'ReadRowNames',false);           
WMt  = readtable(['ds003653' filesep 'WhiteMatter_entropy.csv'],'ReadRowNames',false);           
CSFt = readtable(['ds003653' filesep 'CSF_entropy.csv'],'ReadRowNames',false);   

figure('Name','Tissue Entropy'); 
subplot(3,4,1);
[GMd_est, CId_GM]   = rst_data_plot(GMd{:,:}, 'estimator','trimmed mean','newfig','sub');
title('Gray Matter discovery set','Fontsize',12); ylabel('GM Entropy'); subplot(3,4,5);
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
title('Gray Matter test set','Fontsize',12); ylabel('GM volumes'); subplot(3,4,7);
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

% For vessels, we computed the % of voxels being not gray (<0.1) or gray
% (>0.9), not white (<0.1) or white (>0.9), and not csf (<0.1) or csf (>0.9). 
% We summarize this here by using the ratio, if the ratio is bigger than 1
% it indicates more voxels seen as not from that tissue than from that tissue
% % and conversely, the higher that ratio the better.

GMd       = readtable(['NRU_dataset' filesep 'GrayMatter_vessels.csv'],'ReadRowNames',false);  % High probability of gray matter in vessels
WMd       = readtable(['NRU_dataset' filesep 'WhiteMatter_vessels.csv'],'ReadRowNames',false);           
CSFd      = readtable(['NRU_dataset' filesep 'CSF_vessels.csv'],'ReadRowNames',false);   
notGMd    = readtable(['NRU_dataset' filesep 'Not_GrayMatter_vessels.csv'],'ReadRowNames',false);  % Low probability of gray matter in vessels
notWMd    = readtable(['NRU_dataset' filesep 'Not_WhiteMatter_vessels.csv'],'ReadRowNames',false);           
notCSFd   = readtable(['NRU_dataset' filesep 'Not_CSF_vessels.csv'],'ReadRowNames',false);   
GMratiod  = [mean(notGMd{:,[1 2]} ./ GMd{:,[1 2]},2) mean(notGMd{:,[3 4]} ./ GMd{:,[3 4]},2)];
WMratiod  = [mean(notWMd{:,[1 2]} ./ WMd{:,[1 2]},2) mean(notWMd{:,[3 4]} ./ WMd{:,[3 4]},2)];
CSFratiod = [mean(notCSFd{:,[1 2]}./CSFd{:,[1 2]},2) mean(notCSFd{:,[3 4]}./CSFd{:,[3 4]},2)];

GMt       = readtable(['ds003653' filesep 'GrayMatter_vessels.csv'],'ReadRowNames',false);           
WMt       = readtable(['ds003653' filesep 'WhiteMatter_vessels.csv'],'ReadRowNames',false);           
CSFt      = readtable(['ds003653' filesep 'CSF_vessels.csv'],'ReadRowNames',false);
notGMt    = readtable(['ds003653' filesep 'Not_GrayMatter_vessels.csv'],'ReadRowNames',false);           
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
title('Gray Matter discovery set','Fontsize',12); 
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
title('Gray Matter validation set','Fontsize',12); 
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
title('Not Gray Matter discovery set','Fontsize',12); 
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
title('Not Gray Matter validation set','Fontsize',12); 
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

% where do we see differences spatially
clear variables
Vessels = spm_read_vols(spm_vol(fullfile(fileparts(pwd),...
    ['code' filesep 'Atlases' filesep 'Vessels' filesep 'rvesselRadius.nii'])));
dataset = {'NRU_dataset','ds003653'};
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
    GMd = readtable([dataset{d} filesep 'GrayMatter_volumes.csv'],'ReadRowNames',false);  % High probability of gray matter in vessels
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

disp("Comparing dartel template TIV with the calculated TIV")
disp(spm_summarise(fullfile('nrudataset', filesep, 'sub-52334', filesep, 'anat', filesep, 'templateT1_nG1_6.nii,1'), 'all', 'litres')+ ...
    spm_summarise(fullfile('nrudataset', filesep, 'sub-52334', filesep, 'anat', filesep, 'templateT1_nG1_6.nii,2'), 'all', 'litres')+ ...
    spm_summarise(fullfile('nrudataset', filesep, 'sub-52334', filesep, 'anat', filesep, 'templateT1_nG1_6.nii,3'), 'all', 'litres'))

T = table(spm_summarise(fullfile('nrudataset', filesep, 'sub-52334', filesep, 'anat', filesep, 'templateT1_nG1_6.nii,4'), 'all', 'litres'), ...
          spm_summarise(fullfile('nrudataset', filesep, 'sub-52334', filesep, 'anat', filesep, 'templateT1_nG2_6.nii,4'), 'all', 'litres'), ...
          spm_summarise(fullfile('nrudataset', filesep, 'sub-52334', filesep, 'anat', filesep, 'templateT12_nG1_6.nii,4'), 'all', 'litres'), ...
          spm_summarise(fullfile('nrudataset', filesep, 'sub-52334', filesep, 'anat', filesep, 'templateT12_nG2_6.nii,4'), 'all', 'litres'), ...
          'RowName',{'Volume (L)'}, 'VariableNames', {'T1_nG1', 'T1_nG2', 'T12_nG1', 'T12_nG2'});
      
disp('Discovery set')
disp('-----------------------')
disp(T)
warning('By adding T2w image the tissue class 4 grew with an avage of %g liters',(T.T12_nG1-T.T1_nG1+T.T12_nG2-T.T1_nG2)/2)

% Validation set
%dartel = spm_read_vols(spm_vol(fullfile('ds003653', filesep, 'sub-718216', filesep, 'anat', filesep, 'templateT1_nG1_6.nii,4'))); num_voxels(2,1) = nnz(dartel); % Count the non-zero voxels
%dartel = spm_read_vols(spm_vol(fullfile('ds003653', filesep, 'sub-718216', filesep, 'anat', filesep, 'templateT1_nG2_6.nii,4'))); num_voxels(2,2) = nnz(dartel); % Count the non-zero voxels
%dartel = spm_read_vols(spm_vol(fullfile('ds003653', filesep, 'sub-718216', filesep, 'anat', filesep, 'templateT12_nG1_6.nii,4')));num_voxels(2,3) = nnz(dartel); % Count the non-zero voxels
%dartel = spm_read_vols(spm_vol(fullfile('ds003653', filesep, 'sub-718216', filesep, 'anat', filesep, 'templateT12_nG2_6.nii,4')));num_voxels(2,4) = nnz(dartel); % Count the non-zero voxels

T = table(spm_summarise(fullfile('ds003653', filesep, 'sub-718216', filesep, 'anat', filesep, 'templateT1_nG1_6.nii,4'), 'all', 'litres'), ...
          spm_summarise(fullfile('ds003653', filesep, 'sub-718216', filesep, 'anat', filesep, 'templateT1_nG2_6.nii,4'), 'all', 'litres'), ...
          spm_summarise(fullfile('ds003653', filesep, 'sub-718216', filesep, 'anat', filesep, 'templateT12_nG1_6.nii,4'), 'all', 'litres'), ...
          spm_summarise(fullfile('ds003653', filesep, 'sub-718216', filesep, 'anat', filesep, 'templateT12_nG2_6.nii,4'), 'all', 'litres'), ...
          'RowName',{'Volume (L)'}, 'VariableNames', {'T1_nG1', 'T1_nG2', 'T12_nG1', 'T12_nG2'});
      
disp('Validation set')
disp('-----------------------')
disp(T)
warning('By adding T2w image the tissue class 4 grew with an avage of %g liters',(T.T12_nG1-T.T1_nG1+T.T12_nG2-T.T1_nG2)/2)
