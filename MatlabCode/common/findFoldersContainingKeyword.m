function oFolders = findFoldersContainingKeyword (folders,keyword)
cnt = 1;
for ii = 1:length(folders)
    if folders(ii).isdir
        thisName = lower(folders(ii).name);
        if ~isempty(strfind(thisName,lower(keyword)))
            oFolders(cnt) = folders(ii);
            cnt = cnt + 1;
%             thisName
        end
    end
end

