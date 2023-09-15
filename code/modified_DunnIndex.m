function dunnIndex = modified_DunnIndex(varargin)

% The Dunn Index <https://en.wikipedia.org/wiki/Dunn_index> measures the
% separability of clusters in a image. Here, we propose a modification for
% probability tissues derived from a single image, namely gray matter,
% white matter and CSF.
%
% From principle DI = min(diff between clusters)/max(diff within clusters)
% With tissue classes, virtually every voxel belong to each cluster, our
% modification thus operates in three steps:
% 1 - separate on tissue from the others, for instance GM > (WM+CSF)
% 2 - get the max difference from the tissue (i.e. the range of values)
% 3 - get the min difference between the tissue of interest and the sum of
%     the others, for instance min(GM-(WM+CSF)>0)
%
% FORMAT: dunnIndex = modified_DunnIndex(GM,WM,CSF)
% INPUT   GM, WM, CSF are the tissue images eiher as name, header structure
%         (see spm_vol.m) or matrices
% OUTPUT: dunnIndex is a vector of the modified Dunn Index for each tisue
%         class against the others
%
% Cyril Pernet & Marc Cummings - 2023

%% deal with inputs
if nargin < 3
    error('3 input arguments in expected')
else
    GM  = varargin{1};
    WM  = varargin{2};
    CSF = varargin{3};
end

if ischar(GM);  GM  = spm_vol(GM);  end
if ischar(WM);  WM  = spm_vol(WM);  end
if ischar(CSF); CSF = spm_vol(CSF); end

if isstruct(GM);  GM  = spm_read_vols(GM);  end
if isstruct(WM);  WM  = spm_read_vols(WM);  end
if isstruct(CSF); CSF = spm_read_vols(CSF); end

if any(size(GM) ~= size(WM)) || ...
        any(size(GM) ~= size(CSF))
    error('dimension error, data in have different sizes')
end

tissues          = NaN(size(GM,1), size(GM,2), size(GM,3), 3);
tissues(:,:,:,1) = GM;
tissues(:,:,:,2) = WM;
tissues(:,:,:,3) = CSF;
clear GM WM CSF

figopt = 0; % internal checks turning plotting on to vizualize 
            % what is happening: set to another value than 0
if nargin == 4
    figopt = varargin{4};
end

if figopt ~= 0
    C(1,:) = [1 0 0];
    C(2,:) = [0 1 0];
    C(3,:) = [0 0 1];
    C(4,:) = [1 1 0];
    figure; subplot(2,2,1); 
    for t=1:3
        tmp = squeeze(tissues(:,:,:,t));
        histogram(tmp(:),100,'FaceColor',C(t,:)); hold on
    end
    grid on; title('tissue classes')
    axis([0 1 0 50000]); legend({'GM','WM','CSF'})
end

%% compute
dunnIndex = NaN(1,3);
for tissue_class = 1:3
    OthersIdx = find(~ismember(1:3, tissue_class));
    Others    = tissues(:,:,:,OthersIdx(1)) + tissues(:,:,:,OthersIdx(2));
    Mask      = logical(tissues(:,:,:,tissue_class) > Others);
    
    if figopt ~= 0
        figure; 
        for z=1:size(GM,3)
            imagesc(Mask(:,:,z));
            if tissue_class == 1
                title('GM Mask slices',z);
            elseif tissue_class == 2
                title('WM Mask slices',z);
            else
                title('CSF Mask slices',z);
            end
            pause(0.1);
        end
        close(gcf)
    end

    % diff tissue to others only for voxels of prob(tissue) > prob(others) 
    tmp                 = squeeze(tissues(:,:,:,tissue_class));
    tmp                 = tmp(Mask) - Others(Mask); 
    diff_between        = min(tmp);

    if figopt ~= 0
        subplot(2,2,tissue_class+1)
        histogram(tmp(:),100,'FaceColor',C(4,:)); hold on
    end
    
    tmp           = squeeze(tissues(:,:,:,tissue_class));
    tmp           = tmp(Mask);   
    tmp           = unique(tmp); 
    tmp(tmp==0)   = [];
    max_intra     = range(tmp);

    if ~isempty(diff_between) || ~isempty(max_intra)
        dunnIndex(tissue_class) = diff_between/max_intra;
    end

    if figopt ~= 0
        histogram(tmp(:),100,'FaceColor',C(tissue_class,:)); hold on
        grid on; axis([0 1 0 50000]); 
        title('Tissue Distributions: DI',dunnIndex(tissue_class))
        if tissue_class == 1
            legend({'Between','GM'});
        elseif tissue_class == 2
            legend({'Between','WM'});
        else
            legend({'Between','CSF'});
        end
    end
end



