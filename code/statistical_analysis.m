% statistical analysis of multispectral segmentation 
% variables are appended with d for discovery (NRU data N=)
% or with a t for test (ds003653).

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
title('TIV test set','Fontsize',12); ylabel('volumes'); 
xlabel('T1-1G T1-2G T12-1G T12-2G'); subplot(2,2,4)
Data = [TIVt(:,1)-TIVt(:,2), TIVt(:,3)-TIVt(:,4),...
    TIVt(:,1)-TIVt(:,3),TIVt(:,2)-TIVt(:,4)];
[h,CI,p] = rst_1ttest(Data,'estimator','trimmed mean','newfig','no');
TIVt_diff = rst_trimmean(Data);
ylabel('volume differences','Fontsize',10); title('Trimmed mean Differences','Fontsize',12);
xlabel('1vs.2 3vs.4 1vs.3 2vs.4')
disp('-----');
warning('test set confirms differences observed in the discovery set');
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

% replication set
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

disp('-------------------------')
disp('     Discovery set')
disp(DI_TM_table)
disp(GMd_table)
disp(WMd_table)
disp(CSFd_table)
disp('-------------------------')

% replication set
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
disp('     Replication set')
disp(DI_table)
disp(GMt_table)
disp(WMt_table)
disp(CSFt_table)
disp('-------------------------')

clearvars

% entropy
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

disp('-------------------------')
disp('     Discovery set')
disp(DI_table)
disp(GMd_table)
disp(WMd_table)
disp(CSFd_table)
disp('-------------------------')

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

GMt_table = table(GMt_dCI(:,:)',GMt_p','VariableNames',{'CI','p-value'},'RowName',{'T1nG1 vs T1nG2', 'T12nG1 vs T12nG2', 'T1nG1 vs T12nG1', 'T1nG2 vs T12nG2'});
WMt_table = table(WMt_dCI(:,:)',WMt_p','VariableNames',{'CI','p-value'},'RowName',{'T1nG1 vs T1nG2', 'T12nG1 vs T12nG2', 'T1nG1 vs T12nG1', 'T1nG2 vs T12nG2'});
CSFt_table = table(CSFt_dCI(:,:)',CSFt_p','VariableNames',{'CI','p-value'},'RowName',{'T1nG1 vs T1nG2', 'T12nG1 vs T12nG2', 'T1nG1 vs T12nG1', 'T1nG2 vs T12nG2'});

disp('-------------------------')
disp('     Replication set')
disp(DI_table)
disp(GMt_table)
disp(WMt_table)
disp(CSFt_table)
disp('-------------------------')

clearvars

% -------------------------------------------------------------------------
%% How much of the missing volumes are now vessels? 
% we investigated where missing volumes are located, thus focusing analyses
% on comparing modality changes only (ie T1w vs, T1w and T2w image inputs in 
% the 1Gaussian and in the 2 Gaussians models)
% -------------------------------------------------------------------------

GMd  = readtable(['NRU_dataset' filesep 'GrayMatter_vessels.csv'],'ReadRowNames',false);           
WMd  = readtable(['NRU_dataset' filesep 'WhiteMatter_vessels.csv'],'ReadRowNames',false);           
CSFd = readtable(['NRU_dataset' filesep 'CSF_vessels.csv'],'ReadRowNames',false);           

GMt  = readtable(['ds003653' filesep 'GrayMatter_vessels.csv'],'ReadRowNames',false);           
WMt  = readtable(['ds003653' filesep 'WhiteMatter_vessels.csv'],'ReadRowNames',false);           
CSFt = readtable(['ds003653' filesep 'CSF_vessels.csv'],'ReadRowNames',false);        

% Discovery set
figure('Name','Tissue vessels'); 
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

% replication set
figure(findobj( 'Type', 'Figure', 'Name', 'Tissue vessels' ));
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




% Display the correlation matrices for each condition
% summery_corrd = table(table(d_correlation_matric(:,1,1), d_correlation_matric(:,2,1), d_correlation_matric(:,3,1),'VariableNames',{'GM', 'WM', 'CSF'}), ...
%     table(d_correlation_matric(:,1,2), d_correlation_matric(:,2,2), d_correlation_matric(:,3,2),'VariableNames',{'GM', 'WM', 'CSF'}), ...
%     table(d_correlation_matric(:,1,3), d_correlation_matric(:,2,3), d_correlation_matric(:,3,3),'VariableNames',{'GM', 'WM', 'CSF'}), ...
%     table(d_correlation_matric(:,1,4), d_correlation_matric(:,2,4), d_correlation_matric(:,3,4),'VariableNames',{'GM', 'WM', 'CSF'}), ...
%     'RowNames',{'GM', 'WM', 'CSF'},'VariableNames',{'T1_nG1','T1_nG2','T12_nG1','T12_nG2'});
% summery_corrt = table(table(t_correlation_matric(:,1,1), t_correlation_matric(:,2,1), t_correlation_matric(:,3,1),'VariableNames',{'GM', 'WM', 'CSF'}), ...
%     table(t_correlation_matric(:,1,2), t_correlation_matric(:,2,2), t_correlation_matric(:,3,2),'VariableNames',{'GM', 'WM', 'CSF'}), ...
%     table(t_correlation_matric(:,1,3), t_correlation_matric(:,2,3), t_correlation_matric(:,3,3),'VariableNames',{'GM', 'WM', 'CSF'}), ...
%     table(t_correlation_matric(:,1,4), t_correlation_matric(:,2,4), t_correlation_matric(:,3,4),'VariableNames',{'GM', 'WM', 'CSF'}), ...
%     'RowNames',{'GM', 'WM', 'CSF'},'VariableNames',{'T1_nG1','T1_nG2','T12_nG1','T12_nG2'});
% disp("Discovery set");
% disp(summery_corrd);
% disp("Test set");
% disp(summery_corrt);


% Add x-axis labels
% XL       = get(findall(figure(1),'type','axes'), 'XLim');
% ticLengh = ((XL(2)-XL(1))/4);
% xticks(ticLengh-XL(1) : ticLengh : (ticLengh*4)-XL(1));
% xticklabels(Conditions);
% XL       = get(findall(figure(2),'type','axes'), 'XLim');
% ticLengh = ((XL(2)-XL(1))/4);
% xticks(ticLengh-XL(1) : ticLengh : (ticLengh*4)-XL(1));
% xticklabels(Conditions);
% XL       = get(findall(figure(3),'type','axes'), 'XLim');
% ticLengh = ((XL(2)-XL(1))/4);
% xticks(ticLengh-XL(1) : ticLengh : (ticLengh*4)-XL(1));
% xticklabels(Conditions);
% if(debug)
%     % save plot for volumes trimmed mean and close figure
%     saveas(figure(1), "GM_Volumes_TM.png"); 
%     saveas(figure(2), "WM_Volumes_TM.png");
%     saveas(figure(3), "CSF_Volumes_TM.png");
% end
% close(figure(1));
% close(figure(2));
% close(figure(3));
% 
% TrimmedMeans = [GM_est; WM_est; CSF_est];
% LowerConfs   = [CI_GM(1,:); CI_WM(1,:); CI_CSF(1,:)];
% HigherConfs  = [CI_GM(2,:); CI_WM(2,:); CI_CSF(2,:)];
% 
% T1_nG1  = [LowerConfs(:,1) TrimmedMeans(:,1) HigherConfs(:,1)];
% T1_nG2  = [LowerConfs(:,2) TrimmedMeans(:,2) HigherConfs(:,2)];
% T12_nG1 = [LowerConfs(:,3) TrimmedMeans(:,3) HigherConfs(:,3)];
% T12_nG2 = [LowerConfs(:,4) TrimmedMeans(:,4) HigherConfs(:,4)];
% 
% T = table(T1_nG1,T1_nG2,T12_nG1,T12_nG2,...
%     'RowNames',TissueNames);
% writetable(T,'Volumes_TrimmedMeans.csv','WriteRowNames',true);
% clear GM_est WM_est CSF_est CI_GM CI_WM CI_CSF TrimmedMeans LowerConfs HigherConfs T1_nG1 T1_nG2 T12_nG1 T12_nG2 T XL ticLengh
% 
% % Multi compare between the 4 conditions (T1_nG1 vs T1_nG2, T12_nG1 vs T12_nG2, T1_nG1 vs T12_nG1, T1_nG2 vs T12_nG2)
% [diff_GM,CI_GM,p_GM,alphav_GM,h_GM]      = rst_multicompare(volumes_GM,[1 2; 3 4; 1 3; 2 4], 'estimator', 'trimmed mean','newfig','yes');
% [diff_WM,CI_WM,p_WM,alphav_WM,h_WM]      = rst_multicompare(volumes_WM,[1 2; 3 4; 1 3; 2 4], 'estimator', 'trimmed mean','newfig','yes');
% [diff_CSF,CI_CSF,p_CSF,alphav_CSF,h_CSF] = rst_multicompare(volumes_CSF,[1 2; 3 4; 1 3; 2 4], 'estimator', 'trimmed mean','newfig','yes');
% % Add x-axis labels
% XL       = get(findall(figure(1),'type','axes'), 'XLim');
% ticLengh = ((XL(2)-XL(1))/4);
% xticks(ticLengh : ticLengh : (ticLengh*4));
% xticklabels(Conditions);
% XL       = get(findall(figure(2),'type','axes'), 'XLim');
% ticLengh = ((XL(2)-XL(1))/4);
% xticks(ticLengh : ticLengh : (ticLengh*4));
% xticklabels(Conditions);
% XL       = get(findall(figure(3),'type','axes'), 'XLim');
% ticLengh = ((XL(2)-XL(1))/4);
% xticks(ticLengh : ticLengh : (ticLengh*4));
% xticklabels(Conditions);
% if(debug)
%     % save plot for Multi compare and close figure
%     saveas(figure(1), "MultiComp_GM.png");
%     saveas(figure(2), "MultiComp_WM.png");
%     saveas(figure(3), "MultiComp_CSF.png");
% end
% close(figure(1));
% close(figure(2));
% close(figure(3));
% 
% PairwiseDifferences = [diff_GM; diff_WM; diff_CSF];
% LowerConfs          = [CI_GM(1,:); CI_WM(1,:); CI_CSF(1,:)];
% HigherConfs         = [CI_GM(2,:); CI_WM(2,:); CI_CSF(2,:)];
% PValues             = [p_GM; p_WM; p_CSF];
% AlphaValues         = [alphav_GM; alphav_WM; alphav_CSF];
% Significances       = [h_GM'; h_WM'; h_CSF'];
% 
% T = table(PairwiseDifferences,LowerConfs,HigherConfs, ...
%           PValues,AlphaValues,Significances,...
%           'RowNames',TissueNames);
% writetable(T,'Volumes_multicompare.csv','WriteRowNames',true);
% clear diff_GM diff_WM diff_CSF CI_GM CI_WM CI_CSF p_GM p_WM p_CSF alphav_GM alphav_WM h_GM h_WM h_CSF alphav_CSF PairwiseDifferences LowerConfs HigherConfs PValues AlphaValues Significances T
% 
% % correlations GM vs WM GM vs CSF WM vs CSF
% [rp_T1_nG1,tp_T1_nG1,CI_T1_nG1,pval_T1_nG1,outid_T1_nG1,h_T1_nG1]       = skipped_Spearman(T1_nG1_vol,[1 2; 1 3; 3 2]);
% [rp_T1_nG2,tp_T1_nG2,CI_T1_nG2,pval_T1_nG2,outid_T1_nG2,h_T1_nG2]       = skipped_Spearman(T1_nG2_vol,[1 2; 1 3; 3 2]);
% [rp_T12_nG1,tp_T12_nG1,CI_T12_nG1,pval_T12_nG1,outid_T12_nG1,h_T12_nG1] = skipped_Spearman(T12_nG1_vol,[1 2; 1 3; 3 2]);
% [rp_T12_nG2,tp_T12_nG2,CI_T12_nG2,pval_T12_nG2,outid_T12_nG2,h_T12_nG2] = skipped_Spearman(T12_nG2_vol,[1 2; 1 3; 3 2]);
% 
% SpearmanCorrelations = [rp_T1_nG1'; rp_T1_nG2'; rp_T12_nG1'; rp_T12_nG2'];
% TValues              = [tp_T1_nG1'; tp_T1_nG2'; tp_T12_nG1'; tp_T12_nG2'];
% LowerConfs           = [CI_T1_nG1(1,:); CI_T1_nG2(1,:); CI_T12_nG1(1,:); CI_T12_nG2(1,:)];
% HigherConfs          = [CI_T1_nG1(2,:); CI_T1_nG2(2,:); CI_T12_nG1(2,:); CI_T12_nG2(2,:)];
% PValues              = [pval_T1_nG1'; pval_T1_nG2'; pval_T12_nG1'; pval_T12_nG2'];
% IndexOutliers        = [outid_T1_nG1'; outid_T1_nG2'; outid_T12_nG1'; outid_T12_nG2'];
% Significances        = [h_T1_nG1'; h_T1_nG2'; h_T12_nG1'; h_T12_nG2'];
% 
% T = table(SpearmanCorrelations,TValues,LowerConfs, ...
%           HigherConfs,PValues,IndexOutliers, Significances,...
%           'RowNames',Conditions);
% writetable(T,'Volumes_spearman.csv','WriteRowNames',true);
% clear rp_T1_nG1 tp_T1_nG1 CI_T1_nG1 pval_T1_nG1 outid_T1_nG1 h_T1_nG1 rp_T1_nG2 tp_T1_nG2 CI_T1_nG2 pval_T1_nG2 outid_T1_nG2 h_T1_nG2 rp_T12_nG1 tp_T12_nG1 CI_T12_nG1 pval_T12_nG1 outid_T12_nG1 h_T12_nG1 rp_T12_nG2 tp_T12_nG2 CI_T12_nG2 pval_T12_nG2 outid_T12_nG2 h_T12_nG2 SpearmanCorrelations TValues LowerConfs HigherConfs PValues IndexOutliers Significances T TissueNames ticLengh XL
% 
% % Total IOntracranial Volume (GM+WM+CSF) -- link to dartel template, where
% % does extra-missing volumes go?
% TIV = volumes_GM + volumes_WM + volumes_CSF;
% [T1_nG1_est, CI_T1_nG1]   = rst_data_plot(TIV(:,1), 'estimator','trimmed mean');
% [T1_nG2_est, CI_T1_nG2]   = rst_data_plot(TIV(:,2), 'estimator','trimmed mean','newfig','yes');
% [T12_nG1_est, CI_T12_nG1] = rst_data_plot(TIV(:,3), 'estimator','trimmed mean','newfig','yes');
% [T12_nG2_est, CI_T12_nG2] = rst_data_plot(TIV(:,4), 'estimator','trimmed mean','newfig','yes');
% % Add x-axis labels
% XL       = get(findall(figure(1),'type','axes'), 'XLim');
% ticLengh = ((XL(2))/2);
% xticks(ticLengh);
% xticklabels(Conditions{1});
% XL       = get(findall(figure(2),'type','axes'), 'XLim');
% ticLengh = ((XL(2))/2);
% xticks(ticLengh);
% xticklabels(Conditions{2});
% XL       = get(findall(figure(3),'type','axes'), 'XLim');
% ticLengh = ((XL(2))/2);
% xticks(ticLengh);
% xticklabels(Conditions{3});
% XL       = get(findall(figure(4),'type','axes'), 'XLim');
% ticLengh = ((XL(2))/2);
% xticks(ticLengh);
% xticklabels(Conditions{4});
% if(debug)
%     % save plot for volumes trimmed mean and close figure
%     saveas(figure(1), "TIV_T1_nG1_Volume_TM.png"); 
%     saveas(figure(2), "TIV_T1_nG2_Volume_TM.png");
%     saveas(figure(3), "TIV_T12_nG1_Volume_TM.png"); 
%     saveas(figure(4), "TIV_T12_nG2_Volume_TM.png");
% end
% close(figure(1));
% close(figure(2));
% close(figure(3));
% close(figure(4));
% 
% TrimmedMeans = [T1_nG1_est, T1_nG2_est, T12_nG1_est, T12_nG2_est];
% LowerConfs   = [CI_T1_nG1(1,:), CI_T1_nG2(1,:), CI_T12_nG1(1,:), CI_T12_nG2(1,:)];
% HigherConfs  = [CI_T1_nG1(2,:), CI_T1_nG2(2,:), CI_T12_nG1(2,:), CI_T12_nG2(2,:)];
% 
% T1_nG1  = [LowerConfs(:,1) TrimmedMeans(:,1) HigherConfs(:,1)];
% T1_nG2  = [LowerConfs(:,2) TrimmedMeans(:,2) HigherConfs(:,2)];
% T12_nG1 = [LowerConfs(:,3) TrimmedMeans(:,3) HigherConfs(:,3)];
% T12_nG2 = [LowerConfs(:,4) TrimmedMeans(:,4) HigherConfs(:,4)];
% 
% T = table(T1_nG1,T1_nG2,T12_nG1,T12_nG2);
% writetable(T,'TIV_Volume_TrimmedMeans.csv','WriteRowNames',true);
% clear T1_nG1_est CI_T1_nG1 T1_nG2_est CI_T1_nG2 T12_nG1_est CI_T12_nG1 T12_nG2_est CI_T12_nG2 TrimmedMeans LowerConfs HigherConfs T1_nG1 T1_nG2 T12_nG1 T12_nG2 T volumes_GM volumes_WM volumes_CSF T1_nG1_vol T1_nG2_vol T12_nG1_vol T12_nG2_vol XL ticLengh TIV
% 
% % distrib_vessels
% load('distrib_vesselsT1_nG1.mat');  T1_nG1_distrib_vessels  = distrib_vessels; clear distrib_vessels
% load('distrib_vesselsT1_nG2.mat');  T1_nG2_distrib_vessels  = distrib_vessels; clear distrib_vessels
% load('distrib_vesselsT12_nG1.mat'); T12_nG1_distrib_vessels = distrib_vessels; clear distrib_vessels
% load('distrib_vesselsT12_nG2.mat'); T12_nG2_distrib_vessels = distrib_vessels; clear distrib_vessels
% 
% TM_T1_nG1  = [trimmean(T1_nG1_distrib_vessels(:,:,1), 5)' trimmean(T1_nG1_distrib_vessels(:,:,2), 5)' trimmean(T1_nG1_distrib_vessels(:,:,3), 5)'];
% TM_T1_nG2  = [trimmean(T1_nG2_distrib_vessels(:,:,1), 5)' trimmean(T1_nG2_distrib_vessels(:,:,2), 5)' trimmean(T1_nG2_distrib_vessels(:,:,3), 5)'];
% TM_T12_nG1 = [trimmean(T12_nG1_distrib_vessels(:,:,1), 5)' trimmean(T12_nG1_distrib_vessels(:,:,2), 5)' trimmean(T12_nG1_distrib_vessels(:,:,3), 5)'];
% TM_T12_nG2 = [trimmean(T12_nG2_distrib_vessels(:,:,1), 5)' trimmean(T12_nG2_distrib_vessels(:,:,2), 5)' trimmean(T12_nG2_distrib_vessels(:,:,3), 5)'];
% 
% TM_GM  = [TM_T1_nG1(:,1) TM_T1_nG2(:,1) TM_T12_nG1(:,1) TM_T12_nG2(:,1)];
% TM_WM  = [TM_T1_nG1(:,2) TM_T1_nG2(:,2) TM_T12_nG1(:,2) TM_T12_nG2(:,2)];
% TM_CSF = [TM_T1_nG1(:,3) TM_T1_nG2(:,3) TM_T12_nG1(:,3) TM_T12_nG2(:,3)];
% 
% [GM_est, CI_GM]   = rst_data_plot(TM_GM, 'estimator','trimmed mean');
% [WM_est, CI_WM]   = rst_data_plot(TM_WM, 'estimator','trimmed mean','newfig','yes');
% [CSF_est, CI_CSF] = rst_data_plot(TM_CSF, 'estimator','trimmed mean','newfig','yes');
% 
% % Add x-axis labels
% XL       = get(findall(figure(1),'type','axes'), 'XLim');
% ticLengh = ((XL(2)-XL(1))/4);
% xticks(ticLengh-XL(1) : ticLengh : (ticLengh*4)-XL(1));
% xticklabels(Conditions);
% XL       = get(findall(figure(2),'type','axes'), 'XLim');
% ticLengh = ((XL(2)-XL(1))/4);
% xticks(ticLengh-XL(1) : ticLengh : (ticLengh*4)-XL(1));
% xticklabels(Conditions);
% XL       = get(findall(figure(3),'type','axes'), 'XLim');
% ticLengh = ((XL(2)-XL(1))/4);
% xticks(ticLengh-XL(1) : ticLengh : (ticLengh*4)-XL(1));
% xticklabels(Conditions);
% if(debug)
%     % save plot for volumes trimmed mean and close figure
%     saveas(figure(1), "GM_vessels_TM.png"); 
%     saveas(figure(2), "WM_vessels_TM.png");
%     saveas(figure(3), "CSF_vessels_TM.png");
% end
% close(figure(1));
% close(figure(2));
% close(figure(3));
% 
% TrimmedMeans = [GM_est; WM_est; CSF_est];
% LowerConfs   = [CI_GM(1,:); CI_WM(1,:); CI_CSF(1,:)];
% HigherConfs  = [CI_GM(2,:); CI_WM(2,:); CI_CSF(2,:)];
% 
% T1_nG1  = [LowerConfs(:,1) TrimmedMeans(:,1) HigherConfs(:,1)];
% T1_nG2  = [LowerConfs(:,2) TrimmedMeans(:,2) HigherConfs(:,2)];
% T12_nG1 = [LowerConfs(:,3) TrimmedMeans(:,3) HigherConfs(:,3)];
% T12_nG2 = [LowerConfs(:,4) TrimmedMeans(:,4) HigherConfs(:,4)];
% 
% T = table(T1_nG1,T1_nG2,T12_nG1,T12_nG2,...
%     'RowNames',TissueNames);
% writetable(T,'Vessels_TrimmedMeans.csv','WriteRowNames',true);
% clear TM_CSF TM_WM TM_GM T1_nG1_distrib_vessels T1_nG2_distrib_vessels T12_nG1_distrib_vessels T12_nG2_distrib_vessels GM_est WM_est CSF_est CI_GM CI_WM CI_CSF TrimmedMeans LowerConfs HigherConfs TM_T1_nG1 TM_T1_nG2 TM_T12_nG1 TM_T12_nG2 T1_nG1 T1_nG2 T12_nG1 T12_nG2 T XL ticLengh
% 
% 
% % HD
% load('HD.mat');
% 
% % entropy
% load('entropyT1_nG1.mat');  T1_nG1_entropy  = entropy; clear entropy
% load('entropyT1_nG2.mat');  T1_nG2_entropy  = entropy; clear entropy
% load('entropyT12_nG1.mat'); T12_nG1_entropy = entropy; clear entropy
% load('entropyT12_nG2.mat'); T12_nG2_entropy = entropy; clear entropy
% 
% entropy_GM  = [T1_nG1_entropy(:,1) T1_nG2_entropy(:,1) T12_nG1_entropy(:,1) T12_nG2_entropy(:,1)];
% entropy_WM  = [T1_nG1_entropy(:,2) T1_nG2_entropy(:,2) T12_nG1_entropy(:,2) T12_nG2_entropy(:,2)];
% entropy_CSF = [T1_nG1_entropy(:,3) T1_nG2_entropy(:,3) T12_nG1_entropy(:,3) T12_nG2_entropy(:,3)];
% 
% [GM_est, CI_GM]   = rst_data_plot(entropy_GM, 'estimator','trimmed mean');
% [WM_est, CI_WM]   = rst_data_plot(entropy_WM, 'estimator','trimmed mean','newfig','yes');
% [CSF_est, CI_CSF] = rst_data_plot(entropy_CSF, 'estimator','trimmed mean','newfig','yes');
% if(debug)
%     % save plot for entropy trimmed mean and close figure
%     saveas(figure(1), "GM_entropy_TM.png"); 
%     saveas(figure(2), "WM_entropy_TM.png");
%     saveas(figure(3), "CSF_entropy_TM.png");
% end
% close(figure(1));
% close(figure(2));
% close(figure(3));
% 
% TrimmedMeans = [GM_est; WM_est; CSF_est];
% LowerConfs   = [CI_GM(1,:); CI_WM(1,:); CI_CSF(1,:)];
% HigherConfs  = [CI_GM(2,:); CI_WM(2,:); CI_CSF(2,:)];
% 
% T1_nG1  = [LowerConfs(:,1) TrimmedMeans(:,1) HigherConfs(:,1)];
% T1_nG2  = [LowerConfs(:,2) TrimmedMeans(:,2) HigherConfs(:,2)];
% T12_nG1 = [LowerConfs(:,3) TrimmedMeans(:,3) HigherConfs(:,3)];
% T12_nG2 = [LowerConfs(:,4) TrimmedMeans(:,4) HigherConfs(:,4)];
% 
% T = table(T1_nG1,T1_nG2,T12_nG1,T12_nG2,...
%     'RowNames',RowNames);
% writetable(T,'entropy_TrimmedMeans.csv','WriteRowNames',true);
% clear GM_est WM_est CSF_est CI_GM CI_WM CI_CSF TrimmedMeans LowerConfs HigherConfs T1_nG1 T1_nG2 T12_nG1 T12_nG2 T
% 
% % Multi compare between the 4 conditions (T1_nG1 vs T1_nG2, T12_nG1 vs T12_nG2, T1_nG1 vs T12_nG1, T1_nG2 vs T12_nG2)
% [diff_GM,CI_GM,p_GM,alphav_GM,h_GM]      = rst_multicompare(entropy_GM,[1 2; 3 4; 1 3; 2 4], 'estimator', 'trimmed mean','newfig','yes');
% [diff_WM,CI_WM,p_WM,alphav_WM,h_WM]      = rst_multicompare(entropy_WM,[1 2; 3 4; 1 3; 2 4], 'estimator', 'trimmed mean','newfig','yes');
% [diff_CSF,CI_CSF,p_CSF,alphav_CSF,h_CSF] = rst_multicompare(entropy_CSF,[1 2; 3 4; 1 3; 2 4], 'estimator', 'trimmed mean','newfig','yes');
% if(debug)
%     % save plot for Multi compare and close figure
%     saveas(figure(1), "MultiComp_GM_entropy.png");
%     saveas(figure(2), "MultiComp_WM_entropy.png");
%     saveas(figure(3), "MultiComp_CSF_entropy.png");
% end
% close(figure(1));
% close(figure(2));
% close(figure(3));
% 
% PairwiseDifferences = [diff_GM; diff_WM; diff_CSF];
% LowerConfs          = [CI_GM(1,:); CI_WM(1,:); CI_CSF(1,:)];
% HigherConfs         = [CI_GM(2,:); CI_WM(2,:); CI_CSF(2,:)];
% PValues             = [p_GM,p_WM, p_CSF];
% AlphaValues         = [alphav_GM; alphav_WM; alphav_CSF];
% Significances       = [h_GM; h_WM; h_CSF];
% 
% clear diff_GM diff_WM diff_CSF CI_GM CI_WM CI_CSF PairwiseDifferences LowerConfs HigherConfs PValues AlphaValues Significances T
% 
% 
% % dunnIndex
% load('dunnIndexT1_nG1.mat');  T1_nG1_dunnIndex  = dunnIndexes; clear dunnIndexes
% load('dunnIndexT1_nG2.mat');  T1_nG2_dunnIndex  = dunnIndexes; clear dunnIndexes
% load('dunnIndexT12_nG1.mat'); T12_nG1_dunnIndex = dunnIndexes; clear dunnIndexes
% load('dunnIndexT12_nG2.mat'); T12_nG2_dunnIndex = dunnIndexes; clear dunnIndexes
% 
% dunnIndex_GM  = [T1_nG1_dunnIndex(:,1) T1_nG2_dunnIndex(:,1) T12_nG1_dunnIndex(:,1) T12_nG2_dunnIndex(:,1)];
% dunnIndex_WM  = [T1_nG1_dunnIndex(:,2) T1_nG2_dunnIndex(:,2) T12_nG1_dunnIndex(:,2) T12_nG2_dunnIndex(:,2)];
% dunnIndex_CSF = [T1_nG1_dunnIndex(:,3) T1_nG2_dunnIndex(:,3) T12_nG1_dunnIndex(:,3) T12_nG2_dunnIndex(:,3)];
% 
% [GM_est, CI_GM]   = rst_data_plot(dunnIndex_GM, 'estimator','trimmed mean');
% [WM_est, CI_WM]   = rst_data_plot(dunnIndex_WM, 'estimator','trimmed mean','newfig','yes');
% [CSF_est, CI_CSF] = rst_data_plot(dunnIndex_CSF, 'estimator','trimmed mean','newfig','yes');
% if(debug)
%     % save plot for dunnIndex trimmed mean and close figure
%     saveas(figure(1), "GM_dunnIndex_TM.png"); 
%     saveas(figure(2), "WM_dunnIndex_TM.png");
%     saveas(figure(3), "CSF_dunnIndex_TM.png");
% end
% close(figure(1));
% close(figure(2));
% close(figure(3));
% 
% TrimmedMeans = [GM_est; WM_est; CSF_est];
% LowerConfs   = [CI_GM(1,:); CI_WM(1,:); CI_CSF(1,:)];
% HigherConfs  = [CI_GM(2,:); CI_WM(2,:); CI_CSF(2,:)];
% 
% T1_nG1  = [LowerConfs(:,1) TrimmedMeans(:,1) HigherConfs(:,1)];
% T1_nG2  = [LowerConfs(:,2) TrimmedMeans(:,2) HigherConfs(:,2)];
% T12_nG1 = [LowerConfs(:,3) TrimmedMeans(:,3) HigherConfs(:,3)];
% T12_nG2 = [LowerConfs(:,4) TrimmedMeans(:,4) HigherConfs(:,4)];
% 
% T = table(T1_nG1,T1_nG2,T12_nG1,T12_nG2,...
%     'RowNames',RowNames);
% writetable(T,'dunnIndex_TrimmedMeans.csv','WriteRowNames',true);
% clear GM_est WM_est CSF_est CI_GM CI_WM CI_CSF TrimmedMeans LowerConfs HigherConfs T1_nG1 T1_nG2 T12_nG1 T12_nG2 T
% 
% % Multi compare between the 4 conditions (T1_nG1 vs T1_nG2, T12_nG1 vs T12_nG2, T1_nG1 vs T12_nG1, T1_nG2 vs T12_nG2)
% [diff_GM,CI_GM,p_GM,alphav_GM,h_GM]      = rst_multicompare(dunnIndex_GM,[1 2; 3 4; 1 3; 2 4], 'estimator', 'trimmed mean','newfig','yes');
% [diff_WM,CI_WM,p_WM,alphav_WM,h_WM]      = rst_multicompare(dunnIndex_WM,[1 2; 3 4; 1 3; 2 4], 'estimator', 'trimmed mean','newfig','yes');
% [diff_CSF,CI_CSF,p_CSF,alphav_CSF,h_CSF] = rst_multicompare(dunnIndex_CSF,[1 2; 3 4; 1 3; 2 4], 'estimator', 'trimmed mean','newfig','yes');
% if(debug)
%     % save plot for Multi compare and close figure
%     saveas(figure(1), "MultiComp_GM_dunnIndex.png");
%     saveas(figure(2), "MultiComp_WM_dunnIndex.png");
%     saveas(figure(3), "MultiComp_CSF_dunnIndex.png");
% end
% close(figure(1));
% close(figure(2));
% close(figure(3));
% 
% PairwiseDifferences = [diff_GM; diff_WM; diff_CSF];
% LowerConfs          = [CI_GM(1,:); CI_WM(1,:); CI_CSF(1,:)];
% HigherConfs         = [CI_GM(2,:); CI_WM(2,:); CI_CSF(2,:)];
% PValues             = [p_GM,p_WM, p_CSF];
% AlphaValues         = [alphav_GM; alphav_WM; alphav_CSF];
% Significances       = [h_GM; h_WM; h_CSF];
% 
% clear diff_GM diff_WM diff_CSF CI_GM CI_WM CI_CSF PairwiseDifferences LowerConfs HigherConfs PValues AlphaValues Significances T
% 
% 
% 
% 
% % distrib
% load('distribT1_nG1.mat');  T1_nG1_distrib  = distrib; clear distrib
% load('distribT1_nG2.mat');  T1_nG2_distrib  = distrib; clear distrib
% load('distribT12_nG1.mat'); T12_nG1_distrib = distrib; clear distrib
% load('distribT12_nG2.mat'); T12_nG2_distrib = distrib; clear distrib
% 
% tissue_distrib_GM  = [T1_nG1_distrib(:,1) T1_nG2_distrib(:,1) T12_nG1_distrib(:,1) T12_nG2_distrib(:,1)];
% tissue_distrib_WM  = [T1_nG1_distrib(:,2) T1_nG2_distrib(:,2) T12_nG1_distrib(:,2) T12_nG2_distrib(:,2)];
% tissue_distrib_CSF = [T1_nG1_distrib(:,3) T1_nG2_distrib(:,3) T12_nG1_distrib(:,3) T12_nG2_distrib(:,3)];
% 
% [GM_est, CI_GM]   = rst_data_plot(tissue_distrib_GM, 'estimator','trimmed mean');
% [WM_est, CI_WM]   = rst_data_plot(tissue_distrib_WM, 'estimator','trimmed mean','newfig','yes');
% [CSF_est, CI_CSF] = rst_data_plot(tissue_distrib_CSF, 'estimator','trimmed mean','newfig','yes');
% if(debug)
%     % save plot for tissue distrib trimmed mean and close figure
%     saveas(figure(1), "GM_tissue_distrib_TM.png"); 
%     saveas(figure(2), "WM_tissue_distrib_TM.png");
%     saveas(figure(3), "CSF_tissue_distrib_TM.png");
% end
% close(figure(1));
% close(figure(2));
% close(figure(3));
% 
% TrimmedMeans = [GM_est; WM_est; CSF_est];
% LowerConfs   = [CI_GM(1,:); CI_WM(1,:); CI_CSF(1,:)];
% HigherConfs  = [CI_GM(2,:); CI_WM(2,:); CI_CSF(2,:)];
% 
% T1_nG1  = [LowerConfs(:,1) TrimmedMeans(:,1) HigherConfs(:,1)];
% T1_nG2  = [LowerConfs(:,2) TrimmedMeans(:,2) HigherConfs(:,2)];
% T12_nG1 = [LowerConfs(:,3) TrimmedMeans(:,3) HigherConfs(:,3)];
% T12_nG2 = [LowerConfs(:,4) TrimmedMeans(:,4) HigherConfs(:,4)];
% 
% T = table(T1_nG1,T1_nG2,T12_nG1,T12_nG2,...
%     'RowNames',RowNames);
% writetable(T,'tissue_distrib_TrimmedMeans.csv','WriteRowNames',true);
% clear GM_est WM_est CSF_est CI_GM CI_WM CI_CSF TrimmedMeans LowerConfs HigherConfs T1_nG1 T1_nG2 T12_nG1 T12_nG2 T
% 
% % Multi compare between the 4 conditions (T1_nG1 vs T1_nG2, T12_nG1 vs T12_nG2, T1_nG1 vs T12_nG1, T1_nG2 vs T12_nG2)
% [diff_GM,CI_GM,p_GM,alphav_GM,h_GM] = rst_multicompare(tissue_distrib_GM,[1 2; 3 4; 1 3; 2 4], 'estimator', 'trimmed mean','newfig','yes');
% [diff_WM,CI_WM,p_WM,alphav_WM,h_WM] = rst_multicompare(tissue_distrib_WM,[1 2; 3 4; 1 3; 2 4], 'estimator', 'trimmed mean','newfig','yes');
% [diff_CSF,CI_CSF,p_CSF,alphav_CSF,h_CSF] = rst_multicompare(tissue_distrib_CSF,[1 2; 3 4; 1 3; 2 4], 'estimator', 'trimmed mean','newfig','yes');
% if(debug)
%     % save plot for Multi compare and close figure
%     saveas(figure(1), "MultiComp_GM_tissue_distrib.png");
%     saveas(figure(2), "MultiComp_WM_tissue_distrib.png");
%     saveas(figure(3), "MultiComp_CSF_tissue_distrib.png");
% end
% close(figure(1));
% close(figure(2));
% close(figure(3));
% 
% PairwiseDifferences = [diff_GM; diff_WM; diff_CSF];
% LowerConfs          = [CI_GM(1,:); CI_WM(1,:); CI_CSF(1,:)];
% HigherConfs         = [CI_GM(2,:); CI_WM(2,:); CI_CSF(2,:)];
% PValues             = [p_GM,p_WM, p_CSF];
% AlphaValues         = [alphav_GM; alphav_WM; alphav_CSF];
% Significances       = [h_GM; h_WM; h_CSF];
% 
% clear diff_GM diff_WM diff_CSF CI_GM CI_WM CI_CSF PairwiseDifferences LowerConfs HigherConfs PValues AlphaValues Significances T
% 
% 
% 
% 
% % distrib_nuclei
% load('distrib_nucleiT1_nG1.mat');  T1_nG1_distrib_nuclei  = distrib_nuclei; clear distrib_nuclei
% load('distrib_nucleiT1_nG2.mat');  T1_nG2_distrib_nuclei  = distrib_nuclei; clear distrib_nuclei
% load('distrib_nucleiT12_nG1.mat'); T12_nG1_distrib_nuclei = distrib_nuclei; clear distrib_nuclei
% load('distrib_nucleiT12_nG2.mat'); T12_nG2_distrib_nuclei = distrib_nuclei; clear distrib_nuclei
% 
% 
% %% Compare derived volumes between the 4 segmentations
% load('volumesT1_nG1.mat');  T1_nG1  = volumes; clear volumes
% load('volumesT1_nG2.mat');  T1_nG2  = volumes; clear volumes
% load('volumesT12_nG1.mat'); T12_nG1 = volumes; clear volumes
% load('volumesT12_nG2.mat'); T12_nG2 = volumes; clear volumes
% 
% %% Compare tissue image distributions between the 4 segmentations
% HD
% distrib
% 
% %% Compare tissue images vozxel-wise between the unimodal and multimodal segmentations with 2 Gaussians
% %t-test
