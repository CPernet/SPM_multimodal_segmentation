function clean_up(fileMap, disable)
%   Clean up unzipped nifti files to spare space on clusters
%
%   Inputs:
%       fileMap: struct('ID',{},'path',{},'T1Path',{},'T2Path',{})
%       disable: logical
%
%   Outputs:
%       None
%
%   Example:
%       clean_up(fileMap, disable)

    if (~disable)
        disp('----------------------------------------')
        disp('   Deleting uncompressed nifti images   ')
        disp('----------------------------------------')
        parfor mapIndex = 1:length(fileMap)
            [filepath,name,ext] = fileparts(fileMap(mapIndex).T1Path);
            delete(string(fileMap(mapIndex).T1Path));
            delete(string(fileMap(mapIndex).T2Path));
            c   = dir([filepath filesep 'c*' T1name]);  for d=1:length(c); delete(fullfile(c(d).folder,c(d).name)); end
            rc  = dir([filepath filesep 'rc*' T1name]); for d=1:length(rc); delete(fullfile(rc(d).folder,rc(d).name)); end
            wc  = dir([filepath filesep 'wc*' T1name]); for d=1:length(wc); delete(fullfile(wc(d).folder,wc(d).name)); end
            seg = dir(append(filepath, filesep, '*seg8.mat')); if ~isempty(seg); delete(fullfile(seg.folder,seg.name)); end
            u_rc = dir(append(filepath, filesep, 'u_rc*')); if ~isempty(u_rc); delete(fullfile(u_rc.folder,u_rc.name)); end
        end
    end
end