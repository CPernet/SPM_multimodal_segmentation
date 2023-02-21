function decompress_gzip(subjectID, path)
%   Decompresses a gzip file
%
%   Inputs:
%       subjectID: ID of subject
%       path: folder path 
%
%   Outputs:
%       None
%
%   Exceptions:
%       IOError: If the file does not exist
%       OSError: If the file is not a gzip file
%
%   Example:
%       decompress_gzip(subjectID, path)

    % Check that the file exists
    if ~isempty(dir(append(path, filesep,subjectID, '*.nii')))
        return;
    end

    % Check that the file is a gzip file
    if isempty(dir(append(path, filesep, subjectID, '*.nii.gz')))
        error('IOError: If the file does not exist');
    end

    % Gunzip the file
    gunzip(append(path, filesep, subjectID, '*.nii.gz'));
end

