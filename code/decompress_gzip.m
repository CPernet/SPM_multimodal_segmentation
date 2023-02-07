function decompress_gzip(datadir,options)
%   Decompresses a gzip file
%   This function decompresses a gzip compressed file.
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
%       decompress_gzip(datadir,options)
    root = cd;
    cd(datadir);
    folders = dir; % maps the subjects' folders
    for subject = 3:length(folders)
        if strcmpi(options.modality,'T1')
            cd(fullfile(datadir, folders(subject).name)+'\newT1');
            % Check that the file exists
            if ~isempty(dir('*.nii'))
                continue;
            end
        
            % Check that the file is a gzip file
            if isempty(dir('*.nii.gz'))
                error('OSError: Not a gzip file');
            end
        
            % Gunzip the file
            gunzip('*.nii.gz');
        elseif strcmpi(options.modality,'T2')
            cd(fullfile(datadir, folders(subject).name)+'\T2w');
            % Check that the file exists
            if exist('*.nii.gz', 'file') ~= 2
                error('IOError: File not found');
            end
        
            % Check that the file is a gzip file
            if isempty(dir('*.nii.gz'))
                error('OSError: Not a gzip file');
            end
        
            % Gunzip the file
            gunzip('*.nii.gz');
        end
    end
    cd(root);
end

