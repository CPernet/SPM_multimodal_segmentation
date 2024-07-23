%% Statistical Analysis of multispectral segmentation 
% 
% The analysis proceed such as NRU data (N=259) are the discovery set, 
% while the public OpenNeuro ds003653 data (N=87) are the validation set.
% For each analysis, we compute effect sizes with confidence intervals to
% evaluate if effects can reproduce and test for significance using NHST to
% evaluate if effects can replicate. After running the demographics and
% data check, any section can be run separatelty.

addpath('external')
cd('../results')

% ----------------------------------------
%% demographics and checking data ordering 
% always run this before any other sections
% ----------------------------------------

% read the data
participantsd = readtable(['nrudataset' filesep 'participants.tsv'], ...
    'ReadRowNames',false,'FileType','delimitedtext');   
participantst = readtable(['ds003653' filesep 'participants.tsv'], ...
    'ReadRowNames',false,'FileType','delimitedtext');   

% load the SPM batch as volumes are reported using this ordering, and 
% ensures that the metadata follow that order allowing to track or group subjects 
checkd = load(['nrudataset' filesep 'batch.mat']); 
checkd = checkd.batch; metad = participantsd;
names  = cellfun(@(x) x{1}.spm.spatial.coreg.estwrite.ref{1}, checkd, 'UniformOutput', false);
ID     = cellfun(@(x) x(min(findstr(x,'sub-')):min(findstr(x,'sub-'))+8), names, 'UniformOutput', false);
for sub = 1:length(ID)
    metad(sub,:) = participantsd(find(cellfun(@(x) strcmp(x,ID{sub}),participantsd.participant_id)),:);
end

checkt = load(['ds003653' filesep 'batch.mat']); 
checkt = checkt.batch; metat = participantst;
names  = cellfun(@(x) x{1}.spm.spatial.coreg.estwrite.ref{1}, checkt, 'UniformOutput', false);
ID     = cellfun(@(x) x(min(findstr(x,'sub-')):min(findstr(x,'sub-'))+9), names, 'UniformOutput', false);
for sub = 1:length(ID)
    metat(sub,:) = participantst(find(cellfun(@(x) strcmp(x,ID{sub}),participantst.participant_id)),:);
end
clear participantsd participantst

% now extract info
STUDIES = unique(metad.study);
scanner = strcmp('ContSSRI',metad.study(:))+(strcmp('NeuroPharm1',metad.study(:))+strcmp('Migræne',metad.study(:)))*2;
fprintf('The discovery dataset has %g subjects\n',size(metad,1))
fprintf('With an average Male/Female sex ratio of %g (%g%% males, %g%% females)\n ', ...
    sum(strcmpi(metad.Gender,'M')) / sum(strcmpi(metad.Gender,'F')), ...
    sum(strcmpi(metad.Gender,'M'))/size(metad,1)*100, ...
    sum(strcmpi(metad.Gender,'F'))/size(metad,1)*100)
fprintf('A mean age of %g y.o. (std %g)\n', mean(metad.Age),std(metad.Age))
fprintf('It is made up from %g studies looking at depression\n',size(STUDIES,1))
fprintf('With an average Control/Patient ratio of %g (%g%% controls, %g%% patients)\n ', ...
    sum(strcmpi(metad.Group,'Cont')) / sum(strcmpi(metad.Group,'MDD')), ...
    sum(strcmpi(metad.Group,'Cont'))/size(metad,1)*100, ...
    sum(strcmpi(metad.Group,'MDD'))/size(metad,1)*100)
for s = 1:size(STUDIES,1)
    index = strcmp(STUDIES{s},metad.study(:));
    fprintf('Study %s, N=%g\n',STUDIES{s},sum(index))
    fprintf('Study %s, Mean age=%g [%g %g]\n',STUDIES{s}, ...
        mean(metad.Age(index)),min(metad.Age(index)),max(metad.Age(index)))
    fprintf('Study %s, Male/Female sex ratio %g \n',STUDIES{s}, ...
        sum(strcmpi(metad.Gender(index),'M')) / sum(strcmpi(metad.Gender(index),'F')))
    fprintf('Study %s, Control/Patient ratio of %g \n',STUDIES{s}, ...
        sum(strcmpi(metad.Group(index),'Cont')) / sum(strcmpi(metad.Group(index),'MDD')))    
end

fprintf('The validation dataset has %g subjects\n',size(metat,1))
fprintf('With an average Male/Female sex ratio of %g (%g%% males, %g%% females)\n ', ...
    sum(strcmpi(metat.sex,'M')) / sum(strcmpi(metat.sex,'F')), ...
    sum(strcmpi(metat.sex,'M'))/size(metat,1)*100, ...
    sum(strcmpi(metat.sex,'F'))/size(metat,1)*100)
fprintf('A mean age of %g y.o. (std %g)\n', mean(metat.age),std(metat.age))
fprintf('Min %g Max %g\n', min(metat.age),max(metat.age))
fprintf('With an average Control/Patient ratio of %g (%g%% controls, %g%% patients)\n ', ...
    sum(contains(metat.group,'HC')) / sum(strcmpi(metat.group,'UD')), ...
    sum(contains(metat.group,'HC'))/size(metat,1)*100, ...
    sum(strcmpi(metat.group,'UD'))/size(metat,1)*100)

% ---------------------
%% 1. Volume estimates
% ---------------------

GMd  = readtable(['nrudataset' filesep 'GrayMatter_volumes.csv'],'ReadRowNames',false);           
WMd  = readtable(['nrudataset' filesep 'WhiteMatter_volumes.csv'],'ReadRowNames',false);           
CSFd = readtable(['nrudataset' filesep 'CSF_volumes.csv'],'ReadRowNames',false);           
GMt  = readtable(['ds003653' filesep 'GrayMatter_volumes.csv'],'ReadRowNames',false);           
WMt  = readtable(['ds003653' filesep 'WhiteMatter_volumes.csv'],'ReadRowNames',false);           
CSFt = readtable(['ds003653' filesep 'CSF_volumes.csv'],'ReadRowNames',false);           

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
subplot(6,13,[5 6]); plot(K1{1},'b','LineWidth',2); hold on; plot(K1{2},'--b','LineWidth',2); 
grid on; set(gca,'xticklabel',{[]}); title('GM volumes')
subplot(6,13,[13+5 13+6 26+5 26+6]); gscatter([GMd{:,1};GMd{:,2}].*1000,[GMd{:,3};GMd{:,4}].*1000,gp);
grid on; axis('square'); xlabel('T1w'); ylabel('T1w & T2w'); legend('off'); axis([400 1000 400 1000])
subplot(6,13,[13+7 26+7]); plot(fliplr(K1{3}),'--r','LineWidth',2); hold on; 
grid on; plot(fliplr(K1{4}),'r','LineWidth',2); set(gca,'xticklabel',{[]}); camroll(-90)

subplot(6,13,[8 9]); plot(K2{1},'b','LineWidth',2); hold on; plot(K2{2},'--b','LineWidth',2); 
grid on; set(gca,'xticklabel',{[]}); title('WM volumes')
subplot(6,13,[13+8 13+9 26+8 26+9]); gscatter([WMd{:,1};WMd{:,2}].*1000,[WMd{:,3};WMd{:,4}].*1000,gp);
grid on; axis('square'); xlabel('T1w'); legend('off'); axis([300 650 300 650])
subplot(6,13,[13+10 26+10]); plot(fliplr(K2{3}),'--r','LineWidth',2); hold on; 
grid on; plot(fliplr(K2{4}),'r','LineWidth',2); set(gca,'xticklabel',{[]}); camroll(-90)

subplot(6,13,[11 12]); plot(K3{1},'b','LineWidth',2); hold on; plot(K3{2},'--b','LineWidth',2); 
grid on; set(gca,'xticklabel',{[]}); title('CSF volumes')
subplot(6,13,[13+11 13+12 26+11 26+12]); gscatter([CSFd{:,1};CSFd{:,2}].*1000,[CSFd{:,3};CSFd{:,4}].*1000,gp);
grid on; axis('square'); xlabel('T1w'); legend('off'); axis([100 450 100 450])
subplot(6,13,[13+13 26+13]); plot(fliplr(K3{3}),'--r','LineWidth',2); hold on; 
grid on; plot(fliplr(K3{4}),'r','LineWidth',2); set(gca,'xticklabel',{[]}); camroll(-90)

% validatation
gp = [repmat({'1 Gaussian'},size(GMt,1),1);repmat({'2 Gaussians'},size(GMt,1),1)]; 
subplot(6,13,[39+5 39+6]); plot(K4{1},'b','LineWidth',2); hold on; plot(K4{2},'--b','LineWidth',2); 
grid on; set(gca,'xticklabel',{[]}); 
subplot(6,13,[52+5 52+6 65+5 65+6]); gscatter([GMt{:,1};GMt{:,2}].*1000,[GMt{:,3};GMt{:,4}].*1000,gp);
grid on; axis('square'); xlabel('T1w'); ylabel('T1w & T2w'); legend('off'); axis([400 1000 400 1000])
subplot(6,13,[52+7 65+7]); plot(K1{3},'--r','LineWidth',2); hold on; 
grid on; plot(K1{4},'r','LineWidth',2); set(gca,'xticklabel',{[]}); camroll(-90)

subplot(6,13,[39+8 39+9]); plot(K5{1},'b','LineWidth',2); hold on; plot(K5{2},'--b','LineWidth',2); 
grid on; set(gca,'xticklabel',{[]}); 
subplot(6,13,[52+8 52+9 65+8 65+9]); gscatter([WMt{:,1};WMt{:,2}].*1000,[WMt{:,3};WMt{:,4}].*1000,gp);
grid on; axis('square'); xlabel('T1w'); legend('off'); axis([300 650 300 650])
subplot(6,13,[52+10 65+10]); plot(fliplr(K5{3}),'--r','LineWidth',2); hold on; 
grid on; plot(fliplr(K5{4}),'r','LineWidth',2); set(gca,'xticklabel',{[]}); camroll(-90)

subplot(6,13,[39+11 39+12]); plot(K6{1},'b','LineWidth',2); hold on; plot(K6{2},'--b','LineWidth',2); 
grid on; set(gca,'xticklabel',{[]}); 
subplot(6,13,[52+11 52+12 65+11 65+12]); gscatter([CSFt{:,1};CSFt{:,2}].*1000,[CSFt{:,3};CSFt{:,4}].*1000,gp);
grid on; axis('square'); xlabel('T1w'); legend('off'); axis([100 450 100 450])
subplot(6,13,[52+13 65+13]); plot(fliplr(K6{3}),'--r','LineWidth',2); hold on; 
grid on; plot(fliplr(K6{4}),'r','LineWidth',2); set(gca,'xticklabel',{[]}); camroll(-90)

% table 1  
summary = table([TIVd_CI(1,1) TIVd_est(1) TIVd_CI(2,1); TIVt_CI(1,1) TIVt_est(1) TIVt_CI(2,1)],...
    [TIVd_CI(1,2) TIVd_est(2) TIVd_CI(2,2); TIVt_CI(1,2) TIVt_est(2) TIVt_CI(2,2)],...
    [TIVd_CI(1,3) TIVd_est(3) TIVd_CI(2,3); TIVt_CI(1,3) TIVt_est(3) TIVt_CI(2,3)],...
    [TIVd_CI(1,4) TIVd_est(4) TIVd_CI(2,4); TIVt_CI(1,4) TIVt_est(4) TIVt_CI(2,4)],...
    'RowNames',{'TIV discovery','TIV test'},'VariableNames',{'T1_nG1','T1_nG2','T12_nG1','T12_nG2'});
disp(summary); 

condition = {'T1-1G','T1-2G','T12-1G','T12-2G'};
for test = 1:4
    out = intersect(linspace(TIVd_CI(1,test),TIVd_CI(2,test)),linspace(TIVt_CI(1,test),TIVt_CI(2,test)));
    if isempty(out)
        fprintf('non overlap of TIV HDI for %s\n',condition{test})
    else
        warning('overlap of TIV HDI for %s',condition{test})
    end
end

summary = table([CId_GM(1,1) GMd_est(1) CId_GM(2,1); CIt_GM(1,1) GMt_est(1) CIt_GM(2,1)],...
    [CId_GM(1,2) GMd_est(2) CId_GM(2,2); CIt_GM(1,2) GMt_est(2) CIt_GM(2,2)],...
    [CId_GM(1,3) GMd_est(3) CId_GM(2,3); CIt_GM(1,3) GMt_est(3) CIt_GM(2,3)],...
    [CId_GM(1,4) GMd_est(4) CId_GM(2,4); CIt_GM(1,4) GMt_est(4) CIt_GM(2,4)],...
    'RowNames',{'GM discovery','GM test'},'VariableNames',{'T1_nG1','T1_nG2','T12_nG1','T12_nG2'});
disp(summary); 

for test = 1:4
    out = intersect(linspace(CId_GM(1,test),CId_GM(2,test)),linspace(CIt_GM(1,test),CIt_GM(2,test)));
    if isempty(out)
        fprintf('non overlap of GM HDI for %s\n',condition{test})
    else
        warning('overlap of TIV GM for %s',condition{test})
    end
end

summary = table([CId_WM(1,1) WMd_est(1) CId_WM(2,1); CIt_WM(1,1) WMt_est(1) CIt_WM(2,1)],...
    [CId_WM(1,2) WMd_est(2) CId_WM(2,2); CIt_WM(1,2) WMt_est(2) CIt_WM(2,2)],...
    [CId_WM(1,3) WMd_est(3) CId_WM(2,3); CIt_WM(1,3) WMt_est(3) CIt_WM(2,3)],...
    [CId_WM(1,4) WMd_est(4) CId_WM(2,4); CIt_WM(1,4) WMt_est(4) CIt_WM(2,4)],...
    'RowNames',{'WM discovery','WM test'},'VariableNames',{'T1_nG1','T1_nG2','T12_nG1','T12_nG2'});
disp(summary); 

for test = 1:4
    out = intersect(linspace(CId_WM(1,test),CId_WM(2,test)),linspace(CIt_WM(1,test),CIt_WM(2,test)));
    if isempty(out)
        fprintf('non overlap of WM HDI for %s\n',condition{test})
    else
        warning('overlap of TIV WM for %s',condition{test})
    end
end

summary = table([CId_CSF(1,1) CSFd_est(1) CId_CSF(2,1); CIt_CSF(1,1) CSFt_est(1) CIt_CSF(2,1)],...
    [CId_CSF(1,2) CSFd_est(2) CId_CSF(2,2); CIt_CSF(1,2) CSFt_est(2) CIt_CSF(2,2)],...
    [CId_CSF(1,3) CSFd_est(3) CId_CSF(2,3); CIt_CSF(1,3) CSFt_est(3) CIt_CSF(2,3)],...
    [CId_CSF(1,4) CSFd_est(4) CId_CSF(2,4); CIt_CSF(1,4) CSFt_est(4) CIt_CSF(2,4)],...
    'RowNames',{'CSF discovery','CSF test'},'VariableNames',{'T1_nG1','T1_nG2','T12_nG1','T12_nG2'});
disp(summary); 

for test = 1:4
    out = intersect(linspace(CId_CSF(1,test),CId_CSF(2,test)),linspace(CIt_CSF(1,test),CIt_CSF(2,test)));
    if isempty(out)
        fprintf('non overlap of CSF HDI for %s\n',condition{test})
    else
        warning('overlap of TIV CSF for %s',condition{test})
    end
end

% ------------------------------------------------------
%% 2. what the effect on total intracranial volume (TIV) 
% -------------------------------------------------------

% in the discovery set test main effects and interaction using a Hotelling
% test (repeated measure ANOVA) and multiple pair differences (alphav is adjusted
% using Hochberg step-up procedure)

result = rst_rep_anova_T2(TIVd,scanner,[2 2],1000,{'modality','n_gaussians'});
disp('-----');
disp('no scanner effect')
disp(result.gp)
warning('significant effect of the nb of gaussians %g ml and of modality %g ml, with no interaction',...
    mean(rst_trimmean(TIVd(:,[2 4])-TIVd(:,[1 3]))),mean(rst_trimmean(TIVd(:,[3 4])-TIVd(:,[1 2]))))
disp(result.repeated_measure)

warning('significant interaction between scanners and segmentation')
warning('Adding a T2w image using 1 Gaussian or 2 Gaussian has a difference of %g for the Prisma and %g for the Prisma-fit',...
    rst_trimmean(TIVd(scanner==1,3)-TIVd(scanner==1,1)) - rst_trimmean(TIVd(scanner==1,4)-TIVd(scanner==1,2)), ...
    rst_trimmean(TIVd(scanner==2,3)-TIVd(scanner==2,1)) - rst_trimmean(TIVd(scanner==2,4)-TIVd(scanner==2,2)))
disp(result.interaction)
disp('-----')

% replication set - test for the same differences found as above using Bonferonni correction
Data1 = [TIVd(:,2)-TIVd(:,1), TIVd(:,4)-TIVd(:,3),...
    TIVd(:,3)-TIVd(:,1),TIVd(:,4)-TIVd(:,2)];
Data2 = [TIVt(:,2)-TIVt(:,1), TIVt(:,4)-TIVt(:,3),...
    TIVt(:,3)-TIVt(:,1),TIVt(:,4)-TIVt(:,2)];
[h,~,p] = rst_1ttest([mean(Data2(:,[3 4]),2)-mean(Data2(:,[1 2]),2) ...
    mean(Data2(:,[2 4]),2)-mean(Data2(:,[1 3]),2)],'estimator','trimmed mean','figure','off');
disp('-----');
warning('validation set confirms differences observed in the discovery set p<%g',max(p));
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

condition = {'T1 G2-G1','T12 G2-G1','G1 T12-T1','G2 T12-T1'};
for test = 1:4
    out = intersect(linspace(CId_diff(1,test),CId_diff(2,test)),linspace(CIt_diff(1,test),CIt_diff(2,test)));
    if isempty(out)
        fprintf('non overlap of volume difference HDI for %s\n',condition{test})
    else
        warning('overlap of volume difference CSF for %s',condition{test})
    end
end


% find subjects with average effect to illustrate 
% tmp = (abs(Data1) - abs(TIVd_diff));
% [~,index1]=sort([tmp(:,1)+tmp(:,2) tmp(:,3)+tmp(:,4)]);
% tmp = (abs(Data2) - abs(TIVt_diff));
% [~,index2]=sort([tmp(:,1)+tmp(:,2) tmp(:,3)+tmp(:,4)]);
% warning('typical subjects for TIV difference due to parametrization')
% metad(index1(1,1),:)
% metat(index2(1,1),:)
% warning('typical subjects for TIV difference due to input')
% metad(index1(1,2),:)
% metat(index2(1,2),:)


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
fprintf('and %g ml in the ''other'' class\n', ...
(mean(Otherd{:,3}-Otherd{:,1}+Otherd{:,4}-Otherd{:,2})+...
    mean(Othert{:,3}-Othert{:,1}+Othert{:,4}-Othert{:,2}))*250)

fprintf('Adding a Gaussian leads to %g ml\n', (mean(TIVd(:,2)-TIVd(:,1)+TIVd(:,4)-TIVd(:,3))+...
    mean(TIVt(:,2)-TIVt(:,1)+TIVt(:,4)-TIVt(:,3)))/4)
fprintf('this is compensated by a change of %g ml in soft tissue\n', ...
    (mean(SoftTissued{:,2}-SoftTissued{:,1}+SoftTissued{:,4}-SoftTissued{:,3})+...
    mean(SoftTissuet{:,2}-SoftTissuet{:,1}+SoftTissuet{:,4}-SoftTissuet{:,3}))*250)
fprintf('this is compensated by a change of %g ml in bone tisue\n', ...
    (mean(Skulld{:,2}-Skulld{:,1}+Skulld{:,4}-Skulld{:,4})+...
    mean(Skullt{:,2}-Skullt{:,1}+Skullt{:,4}-Skullt{:,4}))*250)
fprintf('and %g ml in the ''other'' class\n', ...
    (mean(Otherd{:,2}-Otherd{:,1}+Otherd{:,4}-Otherd{:,4})+...
    mean(Othert{:,2}-Othert{:,1}+Othert{:,4}-Othert{:,4}))*250)


% ------------------------------------------------------------------------
%% 3. what the effect on the different tissue types and their relationships
% ------------------------------------------------------------------------

figure('Name','Tissue volume differences'); 

% Grey Matter
% ------------
subplot(2,3,1);
Data = [GMd{:,2}-GMd{:,1}, GMd{:,4}-GMd{:,3}, GMd{:,3}-GMd{:,1}, GMd{:,4}-GMd{:,2}].*1000;
[GMd_diff, CIGMd_diff]   = rst_data_plot(Data, 'estimator','trimmed mean','newfig','sub');
title('GM differences discovery set','Fontsize',12);
subplot(2,3,4);
Data = [GMt{:,2}-GMt{:,1}, GMt{:,4}-GMt{:,3}, GMt{:,3}-GMt{:,1},GMt{:,4}-GMt{:,2}].*1000;
[GMt_diff, CIGMt_diff]   = rst_data_plot(Data, 'estimator','trimmed mean','newfig','sub');
title('GM differences validation set','Fontsize',12);

result = rst_rep_anova_T2(GMd{:,:},scanner,[2 2],1000,{'modality','n_gaussians'});
warning('significant effect of modality, nb of gaussians AND interaction') 
[GMdmeans, GMdCI] = rst_rep_anova_plot(GMd{:,:},ones(259,1),[2 2],3);
disp(result.repeated_measure)
disp('-----')
[~,~,GMd_p,~,h1] = rst_multicompare(GMd{:,:}.*1000,[3 1; 4 2], 'estimator', 'trimmed mean','newfig','no');
fprintf('GM volumes differences by adding T2: %g for 1 Gaussians p=%g, %g for 2 Gaussians p=%g\n',...
    GMd_diff(3),GMd_p(1),GMd_diff(4),GMd_p(2))
[h1,CIx,GMd_p] = rst_1ttest((GMt{:,3}-GMt{:,1})-(GMt{:,4}-GMt{:,2})*1000,'trimmean');
fprintf('interation effect does not replicate %g\n',GMd_p)
disp('-----')
disp('no group effect')
disp(result.gp)
disp('-----')
warning('Adding a T2w image using 1 Gaussian or 2 Gaussian has a difference of %g for the Prisma and %g for the Prisma-fit',...
    (rst_trimmean(GMd{scanner==1,3}-GMd{scanner==1,1}) - rst_trimmean(GMd{scanner==1,4}-GMd{scanner==1,2}))*1000, ...
    (rst_trimmean(GMd{scanner==2,3}-GMd{scanner==2,1}) - rst_trimmean(GMd{scanner==2,4}-GMd{scanner==2,2}))*1000)
disp(result.interaction)

% White matter
% -----------
figure(findobj( 'Type', 'Figure', 'Name', 'Tissue volume differences' ));
subplot(2,3,2);
Data = [WMd{:,2}-WMd{:,1}, WMd{:,4}-WMd{:,3}, WMd{:,3}-WMd{:,1}, WMd{:,4}-WMd{:,2}].*1000;
[WMd_diff, CIWMd_diff]   = rst_data_plot(Data, 'estimator','trimmed mean','newfig','sub');
title('WM differences discovery set','Fontsize',12);
subplot(2,3,5);
Data = [WMt{:,2}-WMt{:,1}, WMt{:,4}-WMt{:,3},  WMt{:,3}-WMt{:,1},WMt{:,4}-WMt{:,2}].*1000;
[WMt_diff, CIWMt_diff]   = rst_data_plot(Data, 'estimator','trimmed mean','newfig','sub');
title('WM differences validation set','Fontsize',12);

result = rst_rep_anova_T2(WMd{:,:},scanner,[2 2],1000,{'modality','n_gaussians'});
warning('significant effect of modality, nb of gaussians AND interaction') 
[WMdmeans, WMdCI] = rst_rep_anova_plot(WMd{:,:},ones(259,1),[2 2],3);
disp(result.repeated_measure)
disp('-----')
[~,~,WMd_p,~,h2] = rst_multicompare(WMd{:,:}.*1000,[3 1; 4 2], 'estimator', 'trimmed mean','newfig','no');
fprintf('WM volumes differences by adding T2: %g for 1 Gaussians p=%g, %g for 2 Gaussians p=%g\n', ...
    WMd_diff(3),WMd_p(1),WMd_diff(4),WMd_p(2))
[h2,CIx,WMd_p] = rst_1ttest((WMt{:,3}-WMt{:,1})-(WMt{:,4}-WMt{:,2})*1000,'trimmean');
fprintf('interation effect does not replicate %g\n',WMd_p)
disp('-----')
disp('no group effect')
disp(result.gp)
disp('-----')
warning('Adding a T2w image using 1 Gaussian or 2 Gaussian has a difference of %g for the Prisma and %g for the Prisma-fit',...
    (rst_trimmean(WMd{scanner==1,3}-WMd{scanner==1,1}) - rst_trimmean(WMd{scanner==1,4}-WMd{scanner==1,2}))*1000, ...
    (rst_trimmean(WMd{scanner==2,3}-WMd{scanner==2,1}) - rst_trimmean(WMd{scanner==2,4}-WMd{scanner==2,2}))*1000)
disp(result.interaction)

% CSF
% ---
figure(findobj( 'Type', 'Figure', 'Name', 'Tissue volume differences' ));
subplot(2,3,3);
Data = [CSFd{:,2}-CSFd{:,1}, CSFd{:,4}-CSFd{:,3}, CSFd{:,3}-CSFd{:,1}, CSFd{:,4}-CSFd{:,2}].*1000;
[CSFd_diff, CICSFd_diff]   = rst_data_plot(Data, 'estimator','trimmed mean','newfig','sub');
title('CSF differences discovery set','Fontsize',12);
subplot(2,3,6);
Data = [CSFt{:,2}-CSFt{:,1}, CSFt{:,4}-CSFt{:,3}, CSFt{:,3}-CSFt{:,1},CSFt{:,4}-CSFt{:,2}].*1000;
[CSFt_diff, CICSFt_diff]   = rst_data_plot(Data, 'estimator','trimmed mean','newfig','sub');
title('CSF differences validation set','Fontsize',12);

result = rst_rep_anova_T2(CSFd{:,:},scanner,[2 2],1000,{'modality','n_gaussians'});
warning('significant effect of modality, nb of gaussians AND interaction') 
[CSFdmeans, CSFdCI] = rst_rep_anova_plot(CSFd{:,:},ones(259,1),[2 2],3);
disp(result.repeated_measure)
disp('-----')
[~,~,CSFd_p,~,h3] = rst_multicompare(CSFd{:,:}.*1000,[3 1; 4 2], 'estimator', 'trimmed mean','newfig','no');
fprintf('CSF volumes differences by adding T2: %g for 1 Gaussians p=%g, %g for 2 Gaussians p=%g\n', ...
    CSFd_diff(3),CSFd_p(1),CSFd_diff(4),CSFd_p(2))
[h3,CIx,CSFd_p] = rst_1ttest((WMt{:,3}-WMt{:,1})-(WMt{:,4}-WMt{:,2})*1000,'trimmean');
fprintf('interation effect does not replicate %g\n',CSFd_p)
disp('-----')
disp('no group effect')
disp(result.gp)
disp('-----')
warning('Adding a T2w image has a difference of %g for the Prisma and %g for the Prisma-fit',...
    (rst_trimmean(((CSFd{scanner==1,3}+CSFd{scanner==1,4})./2)-((CSFd{scanner==1,1}+CSFd{scanner==1,2})./2)))*1000, ...
    (rst_trimmean(((CSFd{scanner==2,3}+CSFd{scanner==2,4})./2)-((CSFd{scanner==2,1}+CSFd{scanner==2,2})./2)))*1000)
warning('Adding a Gaussian has a difference of %g for the Prisma and %g for the Prisma-fit',...
    (rst_trimmean(((CSFd{scanner==1,2}+CSFd{scanner==1,4})./2)-((CSFd{scanner==1,1}+CSFd{scanner==1,3})./2)))*1000, ...
    (rst_trimmean(((CSFd{scanner==2,2}+CSFd{scanner==2,4})./2)-((CSFd{scanner==2,1}+CSFd{scanner==2,3})./2)))*1000)
disp(result.interaction)

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

% Multivariate analysis of volumes
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

summary1 = table(mean([gp_11 gp_21 gp_12 gp_22]',2)*100, mean([gp_31 gp_41 gp_32 gp_42]',2)*100', ...
    mean([gp_51 gp_61 gp_52 gp_62]',2)*100, mean([gp_71 gp_81 gp_72 gp_82]',2)*100, ...
    'RowNames',{'CSF+ 1 Gaussian','CSF- 1 Gaussian', 'CSF+ 2 Gaussians', 'CSF- 2 Gaussians'},...
    'VariableNames',{'GM+WM+','GM+WM-','GM-WM+','GM-WM-'});
if single(sum(summary1{[1 2],:}(:)))~= 100 || single(sum(summary1{[3 4],:}(:))) ~= 100
    error('summary percentage does not add up')
end

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

summary2 = table(mean([gp2_11 gp2_21 gp2_12 gp2_22]',2)*100, mean([gp2_31 gp2_41 gp2_32 gp2_42]',2)*100', ...
    mean([gp2_51 gp2_61 gp2_52 gp2_62]',2)*100, mean([gp2_71 gp2_81 gp2_72 gp2_82]',2)*100, ...
    'RowNames',{'CSF+ 1 Gaussian','CSF- 1 Gaussian', 'CSF+ 2 Gaussians', 'CSF- 2 Gaussians'},...
    'VariableNames',{'GM+WM+','GM+WM-','GM-WM+','GM-WM-'});
if single(sum(summary2{[1 2],:}(:)))~= 100 || single(sum(summary1{[3 4],:}(:))) ~= 100
    error('summary percentage does not add up')
end

warning('Discovery dataset'); disp(summary1); 
warning('Validation dataset'); disp(summary2); 

% Now do a hard clustering
% ------------------------

% cluster the discovery dataset with 1 Gaussian parametrization
Datad            = [GMd{:,3}-GMd{:,1},WMd{:,3}-WMd{:,1},CSFd{:,3}-CSFd{:,1}];
[BICS,BESTMODEL] = mbclust(Datad,8);

figure('Name','Guassian Mixture Modelling')
subplot(2,5,1); plotbic(BICS); grid on; box on; axis square; axis([0.5 8.5 3800 4900]);
[class1,uncertainty] = mixclass(Datad,BESTMODEL.pies,BESTMODEL.mus,BESTMODEL.vars);
C   = zeros(length(class1),3);
C(class1==1,1) = 1; % CC
C(class1==2,2) = 1; % green
C(class1==3,3) = 1; % blue
subplot(2,5,2);
scatter3(Datad(:,1),Datad(:,3),Datad(:,2),30,C.*(1-uncertainty)','filled'); 
xlabel('GM'); ylabel('CSF'); zlabel('WM'); axis square; axis([-0.1 0.05 -0.15 0.1 -0.06 0.01])
title(sprintf('Discovery set, 3 clusters\n mean error: %g',mean(uncertainty)))
warning('1 Gaussian parametrization - discovery dataset post-hoc cluster properties')
for c = unique(class1)
    fprintf('class %g: %g Male, %g Female, %g Control, %g Patients, %g Mean age\n',c, ...
        sum(strcmpi(metad(class1'==c,:).Gender,'M')), ...
        sum(strcmpi(metad(class1'==c,:).Gender,'F')), ...
        sum(strcmpi(metad(class1'==c,:).Group,'Cont')), ...
        sum(strcmpi(metad(class1'==c,:).Group,'MDD')), ...
        mean(metad(class1'==c,:).Age));
end

% test model on validation dataset
Datat = [GMt{:,3}-GMt{:,1},WMt{:,3}-WMt{:,1},CSFt{:,3}-CSFt{:,1}];
[class,uncertainty] = mixclass(Datat,BESTMODEL.pies,BESTMODEL.mus,BESTMODEL.vars);
C   = zeros(length(class),3);
C(class==1,1) = 1; % CC
C(class==2,2) = 1; % green
C(class==3,3) = 1; % blue
subplot(2,5,3);
scatter3(Datat(:,1),Datat(:,3),Datat(:,2),30,C.*(1-uncertainty)','filled'); 
xlabel('GM'); ylabel('CSF'); zlabel('WM'); axis square; axis([-0.1 0.05 -0.15 0.1 -0.06 0.01])
title(sprintf('Validation set, same model\n mean error: %g',mean(uncertainty)))

% cluster the validation dataset
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
warning('1 Gaussian parametrization - validation dataset post-hoc cluster properties')
for c = unique(class)
    fprintf('class %g: %g Male, %g Female, %g Control, %g Patients, %g Mean age\n',c, ...
        sum(strcmpi(metat(class'==c,:).sex,'M')), ...
        sum(strcmpi(metat(class'==c,:).sex,'F')), ...
        sum(strcmpi(metat(class'==c,:).group,'HC')), ...
        sum(strcmpi(metat(class'==c,:).group,'UD')), ...
        mean(metat(class'==c,:).age));
end

% cluster the discovery dataset with 2 Gaussian parametrization
Datad            = [GMd{:,4}-GMd{:,2},WMd{:,4}-WMd{:,2},CSFd{:,4}-CSFd{:,2}];
[BICS,BESTMODEL] = mbclust(Datad,8); clustering2Gd = BESTMODEL;

subplot(2,5,6); plotbic(BICS); grid on; box on; axis square; axis([0.5 8.5 3800 4900]);
[class2,uncertainty] = mixclass(Datad,BESTMODEL.pies,BESTMODEL.mus,BESTMODEL.vars);
subplot(2,5,7);
C   = zeros(length(class2),3);
C(class2==3,1) = 1; % CC
C(class2==2,2) = 1; % green
C(class2==1,3) = 1; % blue
scatter3(Datad(:,1),Datad(:,3),Datad(:,2),30,C.*(1-uncertainty)','filled'); 
xlabel('GM'); ylabel('CSF'); zlabel('WM'); axis square; axis([-0.1 0.05 -0.15 0.1 -0.005 0.005])
title(sprintf('Discovery set, 3 clusters\n mean error: %g',mean(uncertainty)))
warning('2 Gaussian parametrization - discovery dataset post-hoc cluster properties')
for c = unique(class2)
    fprintf('class %g: %g Male, %g Female, %g Control, %g Patients, %g Mean age\n',c, ...
        sum(strcmpi(metad(class2'==c,:).Gender,'M')), ...
        sum(strcmpi(metad(class2'==c,:).Gender,'F')), ...
        sum(strcmpi(metad(class2'==c,:).Group,'Cont')), ...
        sum(strcmpi(metad(class2'==c,:).Group,'MDD')), ...
        mean(metad(class2'==c,:).Age));
end

% test model on validation dataset
Datat = [GMt{:,4}-GMt{:,2},WMt{:,4}-WMt{:,2},CSFt{:,4}-CSFt{:,2}];
[class,uncertainty] = mixclass(Datat,BESTMODEL.pies,BESTMODEL.mus,BESTMODEL.vars);
subplot(2,5,8);
C   = zeros(length(class),3);
C(class==3,1) = 1; % CC
C(class==2,2) = 1; % green
C(class==1,3) = 1; % blue
scatter3(Datat(:,1),Datat(:,3),Datat(:,2),30,C.*(1-uncertainty)','filled'); 
xlabel('GM'); ylabel('CSF'); zlabel('WM'); axis square; axis([-0.1 0.05 -0.15 0.1 -0.005 0.005])
title(sprintf('Validation set, same model\n mean error: %g',mean(uncertainty)))

% cluster the validation dataset
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
warning('2 Gaussian parametrization - validation dataset post-hoc cluster properties')
for c = unique(class)
    fprintf('class %g: %g Male, %g Female, %g Control, %g Patients, %g Mean age\n',c, ...
        sum(strcmpi(metat(class'==c,:).sex,'M')), ...
        sum(strcmpi(metat(class'==c,:).sex,'F')), ...
        sum(strcmpi(metat(class'==c,:).group,'HC')), ...
        sum(strcmpi(metat(class'==c,:).group,'UD')), ...
        mean(metat(class'==c,:).age));
end

% we can see that some subjects are different - worth tracking them
% ie subjects 57, 116, 126, 144
outlier_class = intersect(find(class2==3),find(class1==1));
[~,~,up1] = rst_trimci((CSFd{:,3}*1000-CSFd{:,1}*1000)); 
[~,~,up2] = rst_trimci((CSFd{:,4}*1000-CSFd{:,2}*1000));
fprintf('those 4 subjects have CSF changes of %g and %g\nwhile group estimates upper bound are %g and %g\n',...
    mean(CSFd{outlier_class,3}-CSFd{outlier_class,1})*1000, ...
    mean(CSFd{outlier_class,4}-CSFd{outlier_class,2})*1000, up1,up2)

% same main clusters
fprintf('A cluster in the dicovery and validation set are similar when using 2 Gaussians mu=[%g %g %g] vs [%g %g %g]\n',...
    clustering2Gd.mus(:,2)*1000, clustering2Gt.mus(:,1)*1000)
fprintf('sigma= [%g %g %g] vs [%g %g %g]\n', diag(clustering2Gd.vars(:,:,2)*1000), ...
    diag(clustering2Gt.vars(:,:,1)*1000))
fprintf('other validation cluster mu=[%g %g %g] sigma= [%g %g %g]\n', ...
    clustering2Gt.mus(:,2)*1000,diag(clustering2Gt.vars(:,:,2)*1000))

% CCo the grouping per cluster
summary1 = table(mean([gp_11(class2'==2) gp_21(class2'==2) gp_12(class2'==2) gp_22(class2'==2)]',2)*100, ...
    mean([gp_31(class2'==2) gp_41(class2'==2) gp_32(class2'==2) gp_42(class2'==2)]',2)*100', ...
    mean([gp_51(class2'==2) gp_61(class2'==2) gp_52(class2'==2) gp_62(class2'==2)]',2)*100, ...
    mean([gp_71(class2'==2) gp_81(class2'==2) gp_72(class2'==2) gp_82(class2'==2)]',2)*100, ...
    'RowNames',{'CSF+ 1 Gaussian','CSF- 1 Gaussian', 'CSF+ 2 Gaussians', 'CSF- 2 Gaussians'},...
    'VariableNames',{'GM+WM+','GM+WM-','GM-WM+','GM-WM-'});
if sum(summary1{[3 4],:}(:)) ~= 100
    error('summary percentage does not add up')
end

summary2 = table(mean([gp2_11(class'==1) gp2_21(class'==1) gp2_12(class'==1) gp2_22(class'==1)]',2)*100, ...
    mean([gp2_31(class'==1) gp2_41(class'==1) gp2_32(class'==1) gp2_42(class'==1)]',2)*100', ...
    mean([gp2_51(class'==1) gp2_61(class'==1) gp2_52(class'==1) gp2_62(class'==1)]',2)*100, ...
    mean([gp2_71(class'==1) gp2_81(class'==1) gp2_72(class'==1) gp2_82(class'==1)]',2)*100, ...
    'RowNames',{'CSF+ 1 Gaussian','CSF- 1 Gaussian', 'CSF+ 2 Gaussians', 'CSF- 2 Gaussians'},...
    'VariableNames',{'GM+WM+','GM+WM-','GM-WM+','GM-WM-'});
if sum(summary2{[3 4],:}(:)) ~= 100
    error('summary percentage does not add up')
end

summary3 = table(mean([gp2_11(class'==2) gp2_21(class'==2) gp2_12(class'==2) gp2_22(class'==2)]',2)*100, ...
    mean([gp2_31(class'==2) gp2_41(class'==2) gp2_32(class'==2) gp2_42(class'==2)]',2)*100', ...
    mean([gp2_51(class'==2) gp2_61(class'==2) gp2_52(class'==2) gp2_62(class'==2)]',2)*100, ...
    mean([gp2_71(class'==2) gp2_81(class'==2) gp2_72(class'==2) gp2_82(class'==2)]',2)*100, ...
    'RowNames',{'CSF+ 1 Gaussian','CSF- 1 Gaussian', 'CSF+ 2 Gaussians', 'CSF- 2 Gaussians'},...
    'VariableNames',{'GM+WM+','GM+WM-','GM-WM+','GM-WM-'});
if sum(summary2{[3 4],:}(:)) ~= 100
    error('summary percentage does not add up')
end

warning('Discovery dataset main cluster'); disp(summary1); 
warning('Validation dataset reproducible cluster'); disp(summary2); 
warning('Validation dataset the other cluster'); disp(summary3); 


% -------------------------------------------------------------------------
%% 4. similarlity/differences among tissue distributions
% -------------------------------------------------------------------------

HDd   = load(['nrudataset' filesep 'Harrell-Davis-Deciles.mat']);
GMdd  = [squeeze(HDd.HD{3}(:,:,1))-squeeze(HDd.HD{1}(:,:,1)) ...
    squeeze(HDd.HD{4}(:,:,1))-squeeze(HDd.HD{2}(:,:,1))];
WMdd  = [squeeze(HDd.HD{3}(:,:,2))-squeeze(HDd.HD{1}(:,:,2)) ...
    squeeze(HDd.HD{4}(:,:,2))-squeeze(HDd.HD{2}(:,:,2))];
CSFdd = [squeeze(HDd.HD{3}(:,:,3))-squeeze(HDd.HD{1}(:,:,3)) ...
    squeeze(HDd.HD{4}(:,:,3))-squeeze(HDd.HD{2}(:,:,3))];
HDt   = load(['ds003653' filesep 'Harrell-Davis-Deciles.mat']);
GMtd  = [squeeze(HDt.HD{3}(:,:,1))-squeeze(HDt.HD{1}(:,:,1)) ...
    squeeze(HDt.HD{4}(:,:,1))-squeeze(HDt.HD{2}(:,:,1))];
WMtd  = [squeeze(HDt.HD{3}(:,:,2))-squeeze(HDt.HD{1}(:,:,2)) ...
    squeeze(HDt.HD{4}(:,:,2))-squeeze(HDt.HD{2}(:,:,2))];
CSFtd = [squeeze(HDt.HD{3}(:,:,3))-squeeze(HDt.HD{1}(:,:,3)) ...
    squeeze(HDt.HD{4}(:,:,3))-squeeze(HDt.HD{2}(:,:,3))];


for data_type = 1:3

    if data_type == 1
        dd = GMdd; td = GMtd;
        hy = -0.14;
        ax = [0.5 9.5 -0.15 0.1];
    elseif data_type == 2
        dd = WMdd; td = WMtd;
        figure('Name','Shift function WM');
        hy = -0.18;
        ax = [0.5 9.5 -0.20 0.05];
    else
        dd = CSFdd; td = CSFtd;
        figure('Name','Shift function CSF');
        hy = -0.38;
        ax = [0.5 9.5 -0.4 0.2];
    end

    % plot ordering subject per value to have a color gradient (not 100%
    % accurate sum decile but does the job)
    [~,index1] = sort(sum(dd(:,7:8),2));
    [~,index2] = sort(sum(dd(:,16:17),2));
    [~,index3] = sort(sum(td(:,1:9),2));
    [~,index4] = sort(sum(td(:,10:18),2));

    CC          = [zeros(size(dd,1),1) linspace(1,0.1,size(dd,1))' linspace(0.1,1,size(dd,1))'];
    [h,CI,p]    = rst_1ttest(dd(:,1:9),'trimmean','figure','off');
    h           = single(h);
    h(h==0)     = NaN;
    [h2,CI2,p2] = rst_1ttest(dd(:,10:18),'trimmean','figure','off');
    h2          = single(h2);
    h2(h2==0)   = NaN;
    TMd         = rst_trimmean(dd);

    for s = 1: size(dd,1)
        subplot(2,2,1); hold on
        plot(1:9,dd(index1(s),1:9),'Color',CC(s,:));
        subplot(2,2,3); hold on
        plot(1:9,dd(index2(s),10:18),'Color',CC(s,:));
    end

    subplot(2,2,1); grid on; box on; title('1 Gaussian - Discovery dataset'); s=1;
    for i=1:9
        rectangle('Position',[i-0.2,CI(1,i),0.4,CI(2,i)-CI(1,i)],'Curvature',[0.4 0.4],'LineWidth',2,...
            'FaceColor',[0.5 0.5 0.5],'EdgeColor',[0.35 0.35 0.35]);
        plot([i-0.2 i+0.2],[TMd(s) TMd(s)],'LineWidth',3,'Color',[0.35 0.35 0.35]); s=s+1;
    end
    plot(1:9,h*hy,'r*'); axis(ax);

    subplot(2,2,3); grid on; box on; title('2 Gaussians - Discovery dataset')
    for i=1:9
        rectangle('Position',[i-0.2,CI2(1,i),0.4,CI2(2,i)-CI2(1,i)],'Curvature',[0.4 0.4],'LineWidth',2,...
            'FaceColor',[0.5 0.5 0.5],'EdgeColor',[0.35 0.35 0.35]);
        plot([i-0.2 i+0.2],[TMd(s) TMd(s)],'LineWidth',3,'Color',[0.35 0.35 0.35]); s=s+1;
    end
    plot(1:9,h2*hy,'r*'); axis(ax);

    CC          = [zeros(size(td,1),1) linspace(1,0.1,size(td,1))' linspace(0.1,1,size(td,1))'];
    [h,CI,p]    = rst_1ttest(td(:,1:9),'trimmean','figure','off');
    h           = single(h);
    h(h==0)     = NaN;
    [h2,CI2,p2] = rst_1ttest(td(:,10:18),'trimmean','figure','off');
    h2          = single(h2);
    h2(h2==0)   = NaN;
    TMd         = rst_trimmean(td);
    for s = 1: size(td,1)
        subplot(2,2,2); hold on
        plot(1:9,td(index3(s),1:9),'Color',CC(s,:));
        subplot(2,2,4); hold on
        plot(1:9,td(index4(s),10:18),'Color',CC(s,:));
    end

    subplot(2,2,2); grid on; box on; title('1 Gaussian - Validation dataset'); s=1;
    for i=1:9
        rectangle('Position',[i-0.2,CI(1,i),0.4,CI(2,i)-CI(1,i)],'Curvature',[0.4 0.4],'LineWidth',2,...
            'FaceColor',[0.5 0.5 0.5],'EdgeColor',[0.35 0.35 0.35]);
        plot([i-0.2 i+0.2],[TMd(s) TMd(s)],'LineWidth',3,'Color',[0.35 0.35 0.35]); s=s+1;
    end
    plot(1:9,h*hy,'r*'); axis(ax);

    subplot(2,2,4); grid on; box on; title('2 Gaussians - Validation dataset')
    for i=1:9
        rectangle('Position',[i-0.2,CI2(1,i),0.4,CI2(2,i)-CI2(1,i)],'Curvature',[0.4 0.4],'LineWidth',2,...
            'FaceColor',[0.5 0.5 0.5],'EdgeColor',[0.35 0.35 0.35]);
        plot([i-0.2 i+0.2],[TMd(s) TMd(s)],'LineWidth',3,'Color',[0.35 0.35 0.35]); s=s+1;
    end
    plot(1:9,h2*hy,'r*'); axis(ax);

end


% -----------------------------
%% 5. Accuracy testing for ROI
% -----------------------------

% tissue change must be stronger leading to smaller values in voxels with
% arteries located inside the brain too (not just soft tissue and bones)

% For vessels, we computed the % of voxels being not Grey (<0.1) or Grey
% (>0.9), not white (<0.1) or white (>0.9), and not csf (<0.1) or csf (>0.9). 
% We summarize this here by using the ratio, if the ratio is bigger than 1
% it indicates more voxels seen as not from that tissue than from that tissue
% and conversely -- the higher that ratio the better.

GMdv      = readtable(['nrudataset' filesep 'GrayMatter_vessels.csv'],'ReadRowNames',false);  % High probability of Grey matter in vessels
WMdv      = readtable(['nrudataset' filesep 'WhiteMatter_vessels.csv'],'ReadRowNames',false);           
CSFdv     = readtable(['nrudataset' filesep 'CSF_vessels.csv'],'ReadRowNames',false);   
notGMdv   = readtable(['nrudataset' filesep 'Not_GrayMatter_vessels.csv'],'ReadRowNames',false);  % Low probability of Grey matter in vessels
notWMdv   = readtable(['nrudataset' filesep 'Not_WhiteMatter_vessels.csv'],'ReadRowNames',false);           
notCSFdv  = readtable(['nrudataset' filesep 'Not_CSF_vessels.csv'],'ReadRowNames',false);   

GMtv      = readtable(['ds003653' filesep 'GrayMatter_vessels.csv'],'ReadRowNames',false);           
WMtv      = readtable(['ds003653' filesep 'WhiteMatter_vessels.csv'],'ReadRowNames',false);           
CSFtv     = readtable(['ds003653' filesep 'CSF_vessels.csv'],'ReadRowNames',false);
notGMtv   = readtable(['ds003653' filesep 'Not_GrayMatter_vessels.csv'],'ReadRowNames',false);           
notWMtv   = readtable(['ds003653' filesep 'Not_WhiteMatter_vessels.csv'],'ReadRowNames',false);           
notCSFtv  = readtable(['ds003653' filesep 'Not_CSF_vessels.csv'],'ReadRowNames',false);

% test for differences in ratios
GMratiod  = [notGMdv{:,[1 2]} ./ GMdv{:,[1 2]} ,notGMdv{:,[3 4]} ./ GMdv{:,[3 4]}];
GMratiot  = [notGMtv{:,[1 2]} ./ GMtv{:,[1 2]} ,notGMtv{:,[3 4]} ./ GMtv{:,[3 4]}];
result    = rst_rep_anova_T2(GMratiod,scanner,[2 2],1000,{'modality','n_gaussians'});
meansgd   = rst_rep_anova_plot(GMratiod,ones(259,1),[2 2],3);
disp('significant effects')
disp(result.repeated_measure); 

[~,ci,p]  = rst_multicompare(GMratiod,[3 1;4 2], 'estimator', 'trimmed mean','newfig','no');
fprintf('Adding T2w images increases the GM ratio when using 1 Gaussian [%g %g] p=%g but decreases it with 2 Gaussians [%g %g] p=%g\n', ...
    ci(1,1), ci(2,1), p(1), ci(1,2), ci(2,2), p(2));
meansgt   = rst_rep_anova_plot(GMratiot,ones(87,1),[2 2],3);
[h,~,p]   = rst_1ttest((GMratiot(:,3)-GMratiot(:,1))-(GMratiot(:,4)-GMratiot(:,2)),'trimmean');
[~,ci]    = rst_trimmean([GMratiot(:,3)-GMratiot(:,1), GMratiot(:,4)-GMratiot(:,2)]);
fprintf('interation effect replicates p=%g\n although only 1 Gaussian case shows an increase [%g %g] vs [%g %g]', ...
    p, ci(1,1), ci(2,1), ci(1,2), ci(2,2));
disp('-----')
disp('no scanner effect');
disp(result.gp);
disp('but scanner interacts with the metric')
disp(result.interaction)

WMratiod  = [notWMdv{:,[1 2]} ./ WMdv{:,[1 2]} ,notWMdv{:,[3 4]} ./ WMdv{:,[3 4]}];
WMratiot  = [notWMtv{:,[1 2]} ./ WMtv{:,[1 2]} ,notWMtv{:,[3 4]} ./ WMtv{:,[3 4]}];
result    = rst_rep_anova_T2(WMratiod,scanner,[2 2],1000,{'modality','n_gaussians'});
meanswd   = rst_rep_anova_plot(WMratiod,ones(259,1),[2 2],3);
disp('significant effects')
disp(result.repeated_measure); 

[~,ci,p]  = rst_multicompare(WMratiod,[3 1;4 2], 'estimator', 'trimmed mean','newfig','no');
fprintf('Adding T2w images always increases the WM ratio: 1 Gaussian [%g %g] p=%g, 2 Gaussians [%g %g] p=%g\n', ...
    ci(1,1), ci(2,1), p(1), ci(1,2), ci(2,2), p(2));
meanswt   = rst_rep_anova_plot(WMratiot,ones(87,1),[2 2],3);
[h,~,p]   = rst_1ttest((WMratiot(:,3)-WMratiot(:,1))-(WMratiot(:,4)-WMratiot(:,2)),'trimmean');
[~,ci]    = rst_trimmean([WMratiot(:,3)-WMratiot(:,1), WMratiot(:,4)-WMratiot(:,2)]);
fprintf('interation effect replicates p=%g\n although only 1 Gaussian case shows an increase [%g %g] vs [%g %g]', ...
    p, ci(1,1), ci(2,1), ci(1,2), ci(2,2));
disp('-----')
disp('no scanner effect');
disp(result.gp);
disp('but scanner interacts with the metric')
disp(result.interaction)

CSFratiod = [notCSFdv{:,[1 2]} ./ CSFdv{:,[1 2]} ,notCSFdv{:,[3 4]} ./ CSFdv{:,[3 4]}];
CSFratiot = [notCSFtv{:,[1 2]} ./ CSFtv{:,[1 2]} ,notCSFtv{:,[3 4]} ./ CSFtv{:,[3 4]}];
result    = rst_rep_anova_T2(CSFratiod,scanner,[2 2],1000,{'modality','n_gaussians'});
meanscd   = rst_rep_anova_plot(CSFratiod,ones(259,1),[2 2],3);
disp('significant effects')
disp(result.repeated_measure); 

[~,ci,p] = rst_multicompare(CSFratiod,[3 1;4 2], 'estimator', 'trimmed mean','newfig','no');
fprintf('Adding T2w images always increases the CSF ratio: 1 Gaussian [%g %g] p=%g, 2 Gaussians [%g %g] p=%g\n', ...
    ci(1,1), ci(2,1), p(1), ci(1,2), ci(2,2), p(2));
meansct  = rst_rep_anova_plot(CSFratiot,ones(87,1),[2 2],3);
[h,~,p]  = rst_1ttest((CSFratiot(:,3)-CSFratiot(:,1))-(CSFratiot(:,4)-CSFratiot(:,2)),'trimmean');
[~,ci]   = rst_trimmean([CSFratiot(:,3)-CSFratiot(:,1), CSFratiot(:,4)-CSFratiot(:,2)]);
fprintf('interation effect replicates p=%g\n 1 Gaussian case [%g %g], 2 Gausians [%g %g]', ...
    p, ci(1,1), ci(2,1), ci(1,2), ci(2,2));
disp('-----')
disp('no scanner effect');
disp(result.gp);
disp('but scanner interacts with the metric')
disp(result.interaction)

figure('Name','ratio tests'); 
subplot(2,3,1); scatter(GMratiod(:,1),GMratiod(:,3),30,[0 0 1],'filled'); 
hold on; scatter(GMratiot(:,1),GMratiot(:,3),30,[0 1 0],'filled'); 
plot([0 5],[0 5],'k','LineWidth',2); axis([0 5 0 5])
plot([meansgd(1) meansgt(1)],[meansgd(3) meansgt(3)],'r+','LineWidth',3)
xlabel('prob ratio using T1w only'); ylabel(' prob ratio using T1w and T2w');
title('P(GM<.1)/P(GM>.9)'); subtitle('1 Gaussian'); grid on; axis square
subplot(2,3,4); scatter(GMratiod(:,2),GMratiod(:,4),30,[0 0 1],'filled'); 
hold on; scatter(GMratiot(:,2),GMratiot(:,4),30,[0 1 0],'filled');
plot([0 5],[0 5],'k','LineWidth',2); axis([0 5 0 5])
plot([meansgd(2) meansgt(2)],[meansgd(4) meansgt(4)],'r+','LineWidth',3)
xlabel('prob ratio using T1w only'); ylabel(' prob ratio using T1w and T2w');
subtitle('2 Gaussians'); grid on; axis square; 
subplot(2,3,2); scatter(WMratiod(:,1),WMratiod(:,3),30,[0 0 1],'filled'); 
hold on; scatter(WMratiot(:,1),WMratiot(:,3),30,[0 1 0],'filled'); 
plot([0 110],[0 110],'k','LineWidth',2); axis([0 110 0 110])
plot([meanswd(1) meanswt(1)],[meanswd(3) meanswt(3)],'r+','LineWidth',3)
xlabel('prob ratio using T1w only'); ylabel(' prob ratio using T1w and T2w');
title('P(WM<.1)/P(WM>.9)'); subtitle('1 Gaussian'); grid on; axis square
subplot(2,3,5); scatter(WMratiod(:,2),WMratiod(:,4),30,[0 0 1],'filled'); 
hold on; scatter(WMratiot(:,2),WMratiot(:,4),30,[0 1 0],'filled');
plot([0 110],[0 110],'k','LineWidth',2); axis([0 110 0 110])
plot([meanswd(2) meanswt(2)],[meanswd(4) meanswt(4)],'r+','LineWidth',3)
xlabel('prob ratio using T1w only'); ylabel(' prob ratio using T1w and T2w');
subtitle('2 Gaussians'); grid on; axis square;
subplot(2,3,3);scatter(CSFratiod(:,1),CSFratiod(:,3),30,[0 0 1],'filled'); 
hold on; scatter(CSFratiot(:,1),CSFratiot(:,3),30,[0 1 0],'filled'); 
plot([0 30],[0 30],'k','LineWidth',2);axis([0 30 0 30])
plot([meanscd(1) meansct(1)],[meanscd(3) meansct(3)],'r+','LineWidth',3)
xlabel('prob ratio using T1w only'); ylabel(' prob ratio using T1w and T2w');
title('P(CSF<.1)/P(CSF>.9)'); subtitle('1 Gaussian'); grid on; axis square
subplot(2,3,6); scatter(CSFratiod(:,2),CSFratiod(:,4),30,[0 0 1],'filled'); 
hold on; scatter(CSFratiot(:,2),CSFratiot(:,4),30,[0 1 0],'filled');
plot([0 30],[0 30],'k','LineWidth',2); axis([0 30 0 30])
plot([meanscd(2) meansct(2)],[meanscd(4) meansct(4)],'r+','LineWidth',3)
xlabel('prob ratio using T1w only'); ylabel(' prob ratio using T1w and T2w');
subtitle('2 Gaussians'); grid on; axis square; 


% put those ratios back into percentages context, it is the numerator or
% denominator that changes
disp('------------------------------------')
warning('While not capturing the full range of probabilities, threshoding tissues at 0.1 and 0.9, captured %g%% all vessel voxels',...
    mean([mean(mean(notGMdv{:,:}+GMdv{:,:})) ...
    mean(mean(notWMdv{:,:}+WMdv{:,:})) ...
    mean(mean(notCSFdv{:,:}+CSFdv{:,:})) ...
    mean(mean(notGMtv{:,:}+GMtv{:,:})) ...
    mean(mean(notWMtv{:,:}+WMtv{:,:})) ...
    mean(mean(notCSFtv{:,:}+CSFtv{:,:}))]))
fprintf('Across datasets, %g voxels containing large arteries were classified as GM\n', ...
    mean([mean(mean(GMdv{:,:})) mean(mean(GMtv{:,:}))]));
fprintf('%g voxels containing large arteries were classified as WM\n', ...
    mean([mean(mean(WMdv{:,:})) mean(mean(WMtv{:,:}))]));
fprintf('%g voxels containing large arteries were classified as CSF\n', ...
    mean([mean(mean(CSFdv{:,:})) mean(mean(CSFtv{:,:}))]));
fprintf('%g voxels containing large arteries were classified as not GM-WM-CSF\n', ...
    mean([mean(mean(notGMdv{:,:})) mean(mean(notGMtv{:,:})) ...
    mean(mean(notWMdv{:,:}))  mean(mean(notWMtv{:,:}))...
    mean(mean(notCSFdv{:,:})) mean(mean(notCSFtv{:,:}))]));
warning('checking what is driving ratio change - the denominator or the numerator (not tissue)')
warning('1 Gaussian case')
fprintf('In absolute values, adding T2w means %g%% in GM vs %g%% not GM\n',...
    mean(mean(GMdv{:,3}-GMdv{:,1}) + mean(GMtv{:,3}-GMtv{:,1})), ...
    mean(mean(notGMdv{:,3}-notGMdv{:,1}) + mean(notGMtv{:,3}-notGMtv{:,1} )));
fprintf('In absolute values, adding T2w means %g%% in WM vs %g%% not WM\n',...
    mean(mean(WMdv{:,3}-WMdv{:,1}) + mean(WMtv{:,3}-WMtv{:,1})), ...
    mean(mean(notWMdv{:,3}-notWMdv{:,1}) + mean(notWMtv{:,3}-notWMtv{:,1})));
fprintf('In absolute values, adding T2w means %g%% in CSF vs %g%% not CSF\n',...
    mean(mean(CSFdv{:,3}-CSFdv{:,1}) + mean(CSFtv{:,3}-CSFtv{:,1})), ...
    mean(mean(notCSFdv{:,3}-notCSFdv{:,1}) + mean(notCSFtv{:,3}-notCSFtv{:,1})));
warning('2 Gaussians case')
fprintf('In absolute values, adding T2w means %g%% in GM vs %g%% not GM\n',...
    mean(mean(GMdv{:,4}-GMdv{:,2}) + mean(GMtv{:,4}-GMtv{:,2})), ...
    mean(mean(notGMdv{:,4}-notGMdv{:,2}) + mean(notGMtv{:,4}-notGMtv{:,2} )));
fprintf('In absolute values, adding T2w means %g%% in WM vs %g%% not WM\n',...
    mean(mean(WMdv{:,4}-WMdv{:,2}) + mean(WMtv{:,4}-WMtv{:,2})), ...
    mean(mean(notWMdv{:,4}-notWMdv{:,2}) + mean(notWMtv{:,4}-notWMtv{:,2})));
fprintf('In absolute values, adding T2w means %g%% in CSF vs %g%% not CSF\n',...
    mean(mean(CSFdv{:,4}-CSFdv{:,2}) + mean(CSFtv{:,4}-CSFtv{:,2})), ...
    mean(mean(notCSFdv{:,4}-notCSFdv{:,2}) + mean(notCSFtv{:,4}-notCSFtv{:,2})));

% Analysis for nuclei 
% ------------------------------
% The analysis focuses on GM -- for the sake of completness it is also done
% for WM and CFS but it is not relevant here since nuclei are GM. What is
% relevant is the a-priori probability of being GM vs observed. (WM and CSF
% are in the code but not in the article)

colorindex = rst_colour_maps(7); colorindex = colorindex(1:2:7,:);
labels = {'Putamen','Caudate','Nuc. Accumbens','Ext Amygdala',...
    'Globus palludus (ext)','Globus Pallidus (int)', 'Subs. Nigra (compacta)',...
    'Red Nucleus', 'Subs. Nigra (reticulata)', 'Parabranchial Pig. nucl.',...
    'Ventral tegmentum','Ventral Pallidum','Habenular nuclei',...
    'Hypothalamus',' Mammillary Nucleus','Subthalamic nucleus'};

% plot the data as the prob of a tissue as a function of the prob of
% the atlas -- along the way get some estimates to be used for the ratio
% analysis
for dataset=1:2
    if dataset == 1
        tmp = load(['nrudataset' filesep 'distrib_nucleiT1_nG1.mat']); GMnv.T1nG1   = tmp.distrib_nuclei;
        tmp = load(['nrudataset' filesep 'distrib_nucleiT12_nG1.mat']); GMnv.T12nG1 = tmp.distrib_nuclei;
        tmp = load(['nrudataset' filesep 'distrib_nucleiT1_nG2.mat']); GMnv.T1nG2   = tmp.distrib_nuclei;
        tmp = load(['nrudataset' filesep 'distrib_nucleiT12_nG2.mat']); GMnv.T12nG2 = tmp.distrib_nuclei;
    else
        tmp = load(['ds003653' filesep 'distrib_nucleiT1_nG1.mat']); GMnv.T1nG1   = tmp.distrib_nuclei;
        tmp = load(['ds003653' filesep 'distrib_nucleiT12_nG1.mat']); GMnv.T12nG1 = tmp.distrib_nuclei;
        tmp = load(['ds003653' filesep 'distrib_nucleiT1_nG2.mat']); GMnv.T1nG2   = tmp.distrib_nuclei;
        tmp = load(['ds003653' filesep 'distrib_nucleiT12_nG2.mat']); GMnv.T12nG2 = tmp.distrib_nuclei;
    end

    for tissue = 1:3
        if tissue == 1, figure('Name','GM');
        elseif tissue == 2, figure('Name','WM');
        else, figure('Name','CSF'); end

        for nuclei = size(labels,2):-1:1
            tmp = cellfun(@(x) squeeze(size(x,1)), GMnv.T1nG1{nuclei}); 
            tmp(tmp==0)= []; sample(nuclei) = min(tmp); % takes care of the NaN in nuclei 13
            
            % compute the mean value for same number of voxels in each atlas prob.
            tmp  = cellfun(@(x) squeeze(mean(x(randi(sample(nuclei),size(x,1),1),:,:),1)), GMnv.T1nG1{nuclei}, 'UniformOutput', false);
            data(:,:,tissue,1) = cell2mat(cellfun(@(x) x(:,tissue),tmp,'UniformOutput',false));
            [Mp(nuclei,1,:),CIv(nuclei,1,:,:)] = rst_trimmean(squeeze(data(:,:,tissue,1)));
            
            tmp = cellfun(@(x) squeeze(mean(x(randi(sample(nuclei),size(x,1),1),:,:),1)), GMnv.T12nG1{nuclei}, 'UniformOutput', false);
            data(:,:,tissue,2) = cell2mat(cellfun(@(x) x(:,tissue),tmp,'UniformOutput',false));
            [Mp(nuclei,2,:),CIv(nuclei,2,:,:)] = rst_trimmean(squeeze(data(:,:,tissue,2)));
            
            tmp = cellfun(@(x) squeeze(mean(x(randi(sample(nuclei),size(x,1),1),:,:),1)), GMnv.T1nG2{nuclei}, 'UniformOutput', false);
            data(:,:,tissue,3) = cell2mat(cellfun(@(x) x(:,tissue),tmp,'UniformOutput',false));
            [Mp(nuclei,3,:),CIv(nuclei,3,:,:)] = rst_trimmean(squeeze(data(:,:,tissue,3)));
            
            tmp = cellfun(@(x) squeeze(mean(x(randi(sample(nuclei),size(x,1),1),:,:),1)), GMnv.T12nG2{nuclei}, 'UniformOutput', false);
            data(:,:,tissue,4) = cell2mat(cellfun(@(x) x(:,tissue),tmp,'UniformOutput',false));
            [Mp(nuclei,4,:),CIv(nuclei,4,:,:)] = rst_trimmean(squeeze(data(:,:,tissue,4)));

            subplot(4,4,nuclei);
            for cond = 1:4
                plot(1:9,squeeze(Mp(nuclei,cond,:)),'Linewidth',2,'Color',colorindex(cond,:)); hold on
                if nuclei == 13
                    fillhandle = patch([1:6 6:-1:1], [squeeze(CIv(nuclei,cond,1,1:6))',fliplr(squeeze(CIv(nuclei,cond,2,1:6))')], colorindex(cond,:));
                    set(fillhandle,'EdgeColor',colorindex(cond,:),'FaceAlpha',0.2,'EdgeAlpha',0.8);%set edge color
                    fillhandle = patch([8:9 9:-1:8], [squeeze(CIv(nuclei,cond,1,8:9))',fliplr(squeeze(CIv(nuclei,cond,2,8:9))')], colorindex(cond,:));
                    set(fillhandle,'EdgeColor',colorindex(cond,:),'FaceAlpha',0.2,'EdgeAlpha',0.8);%set edge color
                else
                    fillhandle = patch([1:9 9:-1:1], [squeeze(CIv(nuclei,cond,1,:))',fliplr(squeeze(CIv(nuclei,cond,2,:))')], colorindex(cond,:));
                    set(fillhandle,'EdgeColor',colorindex(cond,:),'FaceAlpha',0.2,'EdgeAlpha',0.8);%set edge color
                end
            end

            title(sprintf('%s\n sample=%g voxels',labels{nuclei},sample(nuclei))); 
            grid on; xticks(1:9); axis([0.5 9.5 0 1])
            xticklabels({'.1','.2','.3','.4','.5','.6','.7','.8','.9'});
            if nuclei == 1 || nuclei == 5 || nuclei == 9 || nuclei == 13
                if tissue == 1, ylabel('GM density');
                elseif tissue == 2, ylabel('WM density');
                else, ylabel('CSF density'); end
            end
            if nuclei >= 13
                xlabel('atlas probability')
            end

            if dataset == 1
                if tissue == 1, GMpd{nuclei} = squeeze(data(:,:,tissue,:));
                elseif tissue == 2, WMpd{nuclei} = squeeze(data(:,:,tissue,:));
                else,  CSFpd{nuclei} = squeeze(data(:,:,tissue,:)); end
            else
                if tissue == 1, GMpv{nuclei} = squeeze(data(:,:,tissue,:));
                elseif tissue == 2, WMpv{nuclei} = squeeze(data(:,:,tissue,:));
                else,  CSFpv{nuclei} = squeeze(data(:,:,tissue,:)); end
            end
        end
    end
    clear data
end

% do some ranking 
for nuclei=1:16
    mean_GMpd(nuclei)  = nanmean(GMpd{nuclei}(:));
    mean_WMpd(nuclei)  = nanmean(WMpd{nuclei}(:));
    mean_CSFpd(nuclei) = nanmean(CSFpd{nuclei}(:));
    mean_GMpv(nuclei)  = nanmean(GMpv{nuclei}(:));
    mean_WMpv(nuclei)  = nanmean(WMpv{nuclei}(:));
    mean_CSFpv(nuclei) = nanmean(CSFpv{nuclei}(:));
end
[v,o]  = sort(mean_GMpd,'descend');
[v2,o2]= sort(mean_GMpv,'descend');

% plot the ratios
for dataset=1:2
    figure
    for nuclei=1:16 % o
        if dataset == 1
            ratio = GMpd{o(nuclei)}./(WMpd{o(nuclei)}+CSFpd{o(nuclei)});
        else
            ratio = GMpv{o2(nuclei)}./(WMpv{o2(nuclei)}+CSFpv{o2(nuclei)});
        end
        subplot(4,4,nuclei);
        for cond = 1:4
            [TM,CI] = rst_trimmean(squeeze(ratio(:,:,cond)));
            plot(1:9,TM,'Linewidth',2,'Color',colorindex(cond,:)); hold on
            if nuclei == 13
                fillhandle = patch([1:6 6:-1:1], [CI(1,1:6) fliplr(CI(2,1:6))], colorindex(cond,:));
                set(fillhandle,'EdgeColor',colorindex(cond,:),'FaceAlpha',0.2,'EdgeAlpha',0.8);%set edge color
                fillhandle = patch([8:9 9:-1:8], [CI(1,8:9) fliplr(CI(2,8:9))], colorindex(cond,:));
                set(fillhandle,'EdgeColor',colorindex(cond,:),'FaceAlpha',0.2,'EdgeAlpha',0.8);%set edge color
            else
                fillhandle = patch([1:9 9:-1:1], [CI(1,:) fliplr(CI(2,:))], colorindex(cond,:));
                set(fillhandle,'EdgeColor',colorindex(cond,:),'FaceAlpha',0.2,'EdgeAlpha',0.8);%set edge color
            end
        end
        if dataset == 1
            title(sprintf('%s\n sample=%g voxels',labels{o(nuclei)},sample(o(nuclei))));
        else
            title(sprintf('%s\n sample=%g voxels',labels{o2(nuclei)},sample(o2(nuclei))));
        end
        grid on; xticks(1:9); xticklabels({'.1','.2','.3','.4','.5','.6','.7','.8','.9'});
        if nuclei == 1 || nuclei == 5 || nuclei == 9 || nuclei == 13
            ylabel('Probability ratio');
        end
        if nuclei >= 13
            xlabel('atlas probability')
        end
    end
end
