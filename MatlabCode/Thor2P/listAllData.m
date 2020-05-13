function listAllData(folderName,varargin)
if ~isempty(strfind(folderName,'153465')) | ~isempty(strfind(folderName,'153466'))
    return;
end

if nargin == 2
    testStr = varargin{1};
end

allDirs = dir(folderName);
dirFlags = [allDirs.isdir]';
allDirs = allDirs(dirFlags);
for ii = 1:length(allDirs)
    if strcmp(allDirs(ii).name,'.') | strcmp(allDirs(ii).name,'..')
        continue;
    end
    if allDirs(ii).isdir
        if nargin == 2
            if ~isempty(strfind(sprintf('%s\\%s',folderName,allDirs(ii).name),testStr))
                display(sprintf('%s\\%s',folderName,allDirs(ii).name));
            end
            listAllData(sprintf('%s\\%s',folderName,allDirs(ii).name),testStr);
        else
            display(sprintf('%s\\%s',folderName,allDirs(ii).name));
            listAllData(sprintf('%s\\%s',folderName,allDirs(ii).name));
        end
    end
end
sprintf('\n \n');