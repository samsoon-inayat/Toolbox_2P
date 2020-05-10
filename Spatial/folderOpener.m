function folderOpener(varargin)

p = inputParser;

default_type = 'mainData';
addOptional(p,'type',default_type,@ischar);
parse(p,varargin{:});

ftype = p.Results.type;

[f,cName] = getFolders;
mainDataFolder = f.mainDataFolder; 
tifFolder = f.tifFolder;
mainFolder = f.sDrive;
pmainFolder = f.nsDrive;
pFolder = f.pFolder;

if strcmp(ftype,'mainData')
    winopen(mainDataFolder);
    return;
end

if strcmp(ftype,'processedData')
    winopen(pFolder);
    return;
end
