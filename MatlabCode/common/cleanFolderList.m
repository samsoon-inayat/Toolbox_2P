function folders = cleanFolderList (ofolders,ignoreFolders)
if isempty(ignoreFolders)
    cc = 1;
    for ii = 1:length(ofolders)
        if strcmp(ofolders(ii).name,'.') || strcmp(ofolders(ii).name,'..')
            continue;
        end
        if ~ofolders(ii).isdir
            continue;
        end
        folders(cc) = ofolders(ii);
        cc = cc + 1;
    end
    return;
end
if ischar(ignoreFolders) || iscell(ignoreFolders)
    cc = 1;
    for ii = 1:length(ofolders)
        if strcmp(ofolders(ii).name,'.') || strcmp(ofolders(ii).name,'..')
            continue;
        end
        if ~ofolders(ii).isdir
                continue;
            end
        found = 0;
        for jj = 1:length(ignoreFolders)
            if strcmp(ofolders(ii).name,ignoreFolders{jj})
                found = 1;
                ignoreFolders(jj) = [];
                break;
            end
        end
        if found
            continue
        else
            folders(cc) = ofolders(ii);
            cc = cc + 1;
        end
    end
    return;
end

if isnumeric(ignoreFolders)
    cc = 1;
    for ii = 1:length(ofolders)
        if strcmp(ofolders(ii).name,'.') || strcmp(ofolders(ii).name,'..')
            continue;
        end
        if ~ofolders(ii).isdir
                continue;
        end
         pos = strfind(ofolders(ii).name,'_');
         trialNum = str2num(ofolders(ii).name((pos+1):end));
         if ~isempty(find(ignoreFolders == trialNum))
             continue;
         end
        folders(cc) = ofolders(ii);
        cc = cc + 1;
    end
    return;
end
