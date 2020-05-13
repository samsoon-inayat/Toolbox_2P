function [folders_data] = getFilesData_EV (eDataFolder,stimulus_types,folders_data)
% tfolders = dir(eDataFolder);
% folders = findFoldersContainingKeyword(tfolders,'Evoked');
% 
%     display('Select folder that contains evoked data');
%     [FolderName] = uigetdir(eDataFolder);
% tfolders = dir(eDataFolder);
% folders = tfolders(3:end);
ii = length(folders_data);

% for iii = 1:length(folders)
%     thisName = folders(iii).name;
%     try
%         stimXics = findStimulusType(thisName,stimulus_types);
%         ii = ii + 1;
%     catch
%         continue;
%     end
%     folders_data(ii).name = thisName;
%     folders_data(ii).root_folder = makeName(thisName,eDataFolder);
%     folders_data(ii).stimulus = stimXics;
%     filePath = makeName(folders_data(ii).name,eDataFolder);
%     sFolders = dir(sprintf('%s\\*.tif',filePath));
    sFolders = dir(sprintf('%s\\*.tif',eDataFolder));
    
    hc = 0; lc = 0; nc = 0;
    for jj = 1:length(sFolders)
        thisName = sFolders(jj).name;
        if ~isempty(strfind(thisName,'delta'))
            continue;
        end
        if strcmp(thisName(1:2),'hi') && ~isempty(str2num(thisName(3)))
            hc = hc + 1;
            folders_data(ii).hiTrials{hc} = thisName(1:(end-4));
        end
        if strcmp(thisName(1:2),'lo') && ~isempty(str2num(thisName(3)))
            lc = lc + 1;
            folders_data(ii).loTrials{lc} = thisName(1:(end-4));
        end
        if strcmp(thisName(1:2),'no') && ~isempty(str2num(thisName(3)))
            nc = nc + 1;
            folders_data(ii).noTrials{nc} = thisName(1:(end-4));
        end
    end
% end


% 
% function stimXics = findStimulusType (name,stimTypes)
% % stimConditions = {'post' 'Post'};
% % side = {'Right' 'Left'};
% % for ii = 1:length(stimTypes)
% %     pos = strfind(name,stimTypes{ii});
% %     if ~isempty(pos)
% %         stimXics.stimulus_type = stimTypes{ii};
% %         break;
% %     end
% % end
% % if ~exist('stimXics','var')
% %     n = 0;
% % end
% % stimXics.side = '';
% % for ii = 1:length(side)
% %     pos = strfind(name,side{ii});
% %     if ~isempty(pos)
% %         stimXics.side = side{ii};
% %         break;
% %     end
% % end
% % 
% % % if strcmp(stimXics.side,'')
% % %     stimXics.craniotomy = 'unilateral';
% % % else
% % %     stimXics.craniotomy = 'bilateral';
% % % end
% % 
% % pos = strfind(name,'Trials');
% % stimXics.trials = str2double(name((pos-2):(pos-1)));
% % 
% % pos = strfind(name,'Hz');
% % stimXics.sampling_frequency = str2double(name((pos-3):(pos-1)));
% % 
% % stimXics.strength = '';
% % pos = strfind(name,'mA');
% % if ~isempty(pos)
% %     stimXics.strength = str2double(name((pos-1):(pos-1)));
% %     posu = strfind(name,'_');
% %     poss = posu(find(posu < pos,1,'last'));
% %     stimXics.strengthStr = name((poss+1):(pos+1));
% % end
% % 
% % stimXics.special_conditions = '';
% % for ii = 1:length(stimConditions)
% %     pos = strfind(name,stimConditions{ii});
% %     if ~isempty(pos)
% %         posu = strfind(name,'_');
% %         stimXics.special_conditions = name(posu(end):end);
% %         break;
% %     end
% % end
% 
% posu = strfind(name,'_');
% ps = 1;
% pe = posu(1)-1;
% stimXics{1} = name(ps:pe);
% for ii = 1:length(posu)
%     ps = posu(ii)+1;
%     if ii < length(posu)
%         pe = posu(ii+1)-1;
%     else
%         pe = length(name);
%     end
%     stimXics{ii+1} = name(ps:pe);
% end
% 
% 
% 
% 
