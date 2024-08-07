function fullfilename = decompress_gzip(fullfilename)

%   Decompresses a gzip file
%
%   Inputs: fullfilename (*.nii.gz)
%   Outputs: unzipped file name (*.nii)

% Check that the file exists
if ~exist(fullfilename,'file')
    return
else
    [~,~,ext] = fileparts(fullfilename);
    if strcmpi(ext,'gz')
        % Gunzip the file else nothing and the name is correct
        fullfilename = gunzip(fullfilename);
    end
end

