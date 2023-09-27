function ent = image_entropy(imagein)

% A simple entropy measure of an MRI image
% ent = -sum((voxel value i/N voxels)*log(voxel value i/N voxels))
%
% FORMAT: ent = image_entropy(imagein)
% INPUT   image in is an image eiher as name, header structure
%         (see spm_vol.m) or matrix
% OUTPUT: ent is the image enthropy
%
% Cyril Pernet 

%% deal with input
if ischar(imagein);  imagein  = spm_vol(imagein);  end
if isstruct(imagein);  imagein  = spm_read_vols(imagein);  end

if numel(size(imagein)) ~= 3
    error('3D image in expected')
end

%% compute
mask   = imagein > eps;
values = unique(imagein(mask));
ent    = entropy(imagein(mask));

