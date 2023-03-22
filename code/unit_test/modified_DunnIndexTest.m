function modified_DunnIndexTest

%% get maps from spm
spmdir    = fileparts(which('spm.m'));
template  = fullfile(spmdir,['tpm' filesep 'TPM.nii']);
Vtemplate = spm_vol(template);
if ~exist(fullfile(spmdir,"tpm/TPM_00001.nii"),'file')
    spm_file_split(Vtemplate);
end
GM  = spm_read_vols(spm_vol(fullfile(spmdir,['tpm' filesep 'TPM_00001.nii'])));
WM  = spm_read_vols(spm_vol(fullfile(spmdir,['tpm' filesep 'TPM_00002.nii'])));
CSF = spm_read_vols(spm_vol(fullfile(spmdir,['tpm' filesep 'TPM_00003.nii'])));

%% Part 1 - simply test that inputs and outputs are correct

dunnIndex = zeros(3,3);

try
    % pass in file names
    dunnIndex(1,:) = modified_DunnIndex(fullfile(spmdir,['tpm' filesep 'TPM_00001.nii']),...
        fullfile(spmdir,['tpm' filesep 'TPM_00002.nii']),...
        fullfile(spmdir,['tpm' filesep 'TPM_00003.nii']),1);

    % pass in structures
    dunnIndex(2,:) = modified_DunnIndex(Vtemplate(1),Vtemplate(2),Vtemplate(3));

    % pass in matrices
    dunnIndex(3,:) = modified_DunnIndex(GM,WM,CSF);
catch inerrors
    error('modifled_DunnIndex fail at passing arguments\n%s',inerrors.message)
end

if sum(sum(diff(dunnIndex,1))) ~= 0
    error('for the same data in, as different formats, outputs are different') 
else
    clear dunnIndex
end

%% Part 2 - test against expected values
% GM from 0.1 to 0.6
GM        = double(GM>0.1);
idx       = find(GM);
GMvalues  = randi([1 6],length(idx),1)./10;
GM(idx)   = GMvalues;
GMintra   = range(GM(idx));
% WM from 0.6 to 1
WM        = double(WM>0.1);
idx       = find(WM);
WMvalues  = randi([5 10],length(idx),1)./10;
WM(idx)   = WMvalues;
WMintra   = range(WM(idx));
CSF = zeros(size(CSF));

% diff_between   = 0.1;
% max_within     = 0.5;
expected_value = 0.1/0.5;
dunnIndex = modified_DunnIndex(GM,WM,CSF);
if dunnIndex(1) == dunnIndex(2) && ...
        single(dunnIndex(1)) == single(expected_value)
    disp('----------------------------------------')
    disp('test succesfull returning expected value')
    disp('----------------------------------------')
else
    disp('-------------------------------------------------')
          warning('test failed NOT returning expected value')
    disp('-------------------------------------------------')
end



