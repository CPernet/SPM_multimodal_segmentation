function image_entropyTest

%% get maps from spm
spmdir    = fileparts(which('spm.m'));
template  = fullfile(spmdir,['tpm' filesep 'TPM.nii']);
Vtemplate = spm_vol(template);
if ~exist(fullfile(spmdir,"tpm/TPM_00001.nii"),'file')
    spm_file_split(Vtemplate);
end

%%  simply test that inputs and outputs are correct

entropy = zeros(2,3);

try
    for tissue_class = 1:3
        % pass in file names
        entropy(1,:) = image_entropy(fullfile(spmdir,['tpm' filesep 'TPM_0000' num2str(tissue_class) '.nii']));
    
        % pass in structures
        entropy(2,:) = image_entropy(Vtemplate(tissue_class));
    end
catch inerrors
    error('modifled_DunnIndex fail at passing arguments\n%s',inerrors.message)
end

if sum(sum(diff(entropy,1))) ~= 0
    error('for the same data in, as different formats, outputs are different') 
else
    clear dunnIndex
end


