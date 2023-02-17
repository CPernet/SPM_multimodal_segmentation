function decompress_gzip()
%   Decompresses a gzip file
%
%   Inputs:
%       datadir
%
%   Outputs:
%       None
%
%   Exceptions:
%       IOError: If the file does not exist
%       OSError: If the file is not a gzip file
%
%   Example:
%       decompress_gzip()

    % Check that the file exists
    if ~isempty(dir('*.nii'))
        return;
    end

    % Check that the file is a gzip file
    if isempty(dir('*.nii.gz'))
        error('IOError: If the file does not exist');
    end

    % Gunzip the file
    gunzip('*.nii.gz');
end

