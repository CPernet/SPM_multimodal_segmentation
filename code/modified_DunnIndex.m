function dunnIndex = modified_DunnIndex(subject_path, nvoxels)
    
    %tissues  = NaN(nvoxels,3);
    tissues   = [0.1:0.1:0.6; 0.5:0.1:1; zeros(1 , 6)];
    dunnIndex = NaN(1,3);
    %for tissue_class = 1:3
    %    tmp = dir([subject_path filesep 'wc' num2str(tissue_class) '*.nii']);
    %    tissues(:, tissue_class) = spm_get_data(fullfile(tmp.folder, tmp.name),[x,y,z]');
    %end
    parfor tissue_class = 1:3
        Others = tissues(~ismember(1:3, tissue_class), :);
        Others = Others(1,:) + Others(2,:);
        Mask   = logical(tissues(ismember(1:3, tissue_class), :) > Others);
        tmp    = tissues(tissue_class, Mask) - Others(Mask);
        min_intra_dist = min(tmp(:));
        
        tmp    = unique(tissues(tissue_class,:)); tmp(tmp == 0) = [];
        max_intra_dist = range(tmp);

        dunnIndex(tissue_class) = min_intra_dist/max_intra_dist;
    
    end
end



