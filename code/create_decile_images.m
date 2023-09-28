function HD = create_decile_images(out,outdir,options)
% create a series of 4D images (deciles) for the grey, white and csf tissue
% classes - Harell-Davis estimators are used to find the distribution
% decile and images subsequently thresholded.
%
% INPUTS out is a cell array of the jobs from segmentation
%        outdir where to save the decile images
%        options structure with 'modality', 'NGaussian','Modulate' options
%                (used for naming here)
%
% OUPUTS HD are the Harell-Davis estimators of deciles, size (length(out)-1,9,3)
%        a 4D image is created in the outdir folder
%
% Cyril Pernet - CCBY


% get the mask
spmroot = fileparts(which('spm'));
for p=1:4
    P(p,:) = [spmroot filesep 'tpm' filesep 'TPM.nii,' num2str(p)]; % template 1 to 4
end
V          = spm_vol(P);
M          = single(spm_read_vols(V));
M          = mean(M,4)>0.01;            % mask image (GM,WM,CSF,Meninges)
nvoxels    = sum(M(:));                 % how many in mask voxels
[x,y,z]    = ind2sub(V(1).dim,find(M)); % location of these voxels
V          = V(1);                      % reset to one image for saving later one

HD = NaN(length(out)-1,9,3);

% make images
for tissue_class = 1:3
    
    % read image and get deciles per subjects
    for subject = 1:length(out)-1
        if isfield(out{subject}{1},'tiss')
            file = cell2mat(out{subject}{1}.tiss(tissue_class).wc);     % tissue class image
        else
            file = cell2mat(out{subject}{2}.tiss(tissue_class).wc);
        end
        fprintf('processing subject %g/%g tissue class %g\n', subject, length(out)-1, tissue_class);
        data = spm_get_data(file,[x,y,z]');                       % in mask data
        parfor d=1:9
            xd(d) = get_HD(data',d./10);                          % decile (9 values)
        end
        HD(subject,:,tissue_class) = xd';
        
        img = spm_read_vols(spm_vol(file));                      % threshold image
        if subject == 1
            img1  = (img<=xd(1)) .* M;
            img2  = (img>=xd(1)) .* (img<xd(2)) .*M;
            img3  = (img>=xd(2)) .* (img<xd(3)) .*M;
            img4  = (img>=xd(3)) .* (img<xd(4)) .*M;
            img5  = (img>=xd(4)) .* (img<xd(5)) .*M;
            img6  = (img>=xd(5)) .* (img<xd(6)) .*M;
            img7  = (img>=xd(6)) .* (img<xd(7)) .*M;
            img8  = (img>=xd(7)) .* (img<xd(8)) .*M;
            img9  = (img>=xd(8)) .* (img<xd(9)) .*M;
            img10 = (img>=xd(9)) .* M;
        else
            img1  = img1  + (img<xd(1))  .* M;
            img2  = img2  + (img>=xd(1)) .* (img<xd(2)) .*M;
            img3  = img3  + (img>=xd(2)) .* (img<xd(3)) .*M;
            img4  = img4  + (img>=xd(3)) .* (img<xd(4)) .*M;
            img5  = img5  + (img>=xd(4)) .* (img<xd(5)) .*M;
            img6  = img6  + (img>=xd(5)) .* (img<xd(6)) .*M;
            img7  = img7  + (img>=xd(6)) .* (img<xd(7)) .*M;
            img8  = img8  + (img>=xd(7)) .* (img<xd(8)) .*M;
            img9  = img9  + (img>=xd(8)) .* (img<xd(9)) .*M;
            img10 = img10 + (img>=xd(9)) .* M;
        end
    end
    
    % create the percentage map
    img1 = img1/subject.*100; img2 = img2/subject.*100;
    img3 = img3/subject.*100; img4 = img4/subject.*100;
    img5 = img5/subject.*100; img6 = img6/subject.*100;
    img7 = img7/subject.*100; img8 = img8/subject.*100;
    img9 = img9/subject.*100; img10 = img10/subject.*100;
    
    % save to disk
    for f=1:10
        V.fname = char(fullfile(outdir,filesep,['decile' num2str(f) options.modality '_nG' num2str(options.NGaussian) '.nii']));
        s(f)    = spm_write_vol(V,eval(['img' num2str(f)]));
    end
    
    if tissue_class == 1
        spm_file_merge(s,char(fullfile(outdir,filesep,['GM_deciles_' options.modality '_nG' num2str(options.NGaussian) '.nii'])))
    elseif tissue_class == 2
        spm_file_merge(s,char(fullfile(outdir,filesep,['WM_deciles_' options.modality '_nG' num2str(options.NGaussian) '.nii'])))
    elseif tissue_class == 3
        spm_file_merge(s,char(fullfile(outdir,filesep,['CSF_deciles_' options.modality '_nG' num2str(options.NGaussian) '.nii'])))
    end
    
    for f=1:10
        spm_unlink(s(f).fname);
    end
end

end

function HD = get_HD(x,q)
% that's the Harrell-Davis estimator
n  = length(x);
m1 = (n+1).*q;
m2 = (n+1).*(1-q);
vec= 1:n;
w  = betacdf(vec./n,m1,m2)-betacdf((vec-1)./n,m1,m2);
y  = sort(x);
HD = sum(w(:).*y(:));
end
