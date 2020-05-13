function [folders_data] = getFilesData_SP (eDataFolder,varargin)
tfolders = dir(eDataFolder);
folders = findFilesContainingKeyword(tfolders,'Spon');
if isempty(folders)
    display('Select multiple tif files');
    [FileName,PathName,~] = uigetfile(sprintf('%s\\.tif',eDataFolder),'MultiSelect','on');
    for ii = 1:length(FileName)
        folders(ii).name = FileName{ii};
    end
    eDataFolder = PathName(1:(end-1));
end

folders_data(length(folders)).name = '';
for ii = 1:length(folders)
    thisName = folders(ii).name;
    folders_data(ii).name = thisName(1:(end-4));
    folders_data(ii).root_folder = eDataFolder;
    folders_data(ii).props = findProps(thisName(1:(end-4)));
end


function props = findProps (name)
% 
% pos = strfind(lower(name),'hz');
% if isempty(pos)
% %     display(name);
% %     xx = input('Enter Sampling Frequency = ');
%     props.sampling_frequency = -1;
% else
%     props.sampling_frequency = str2double(name((pos-3):(pos-1)));
% end
% 
% 
% props.number_of_frames = [];
% pos = strfind(lower(name),'fr');
% if length(pos) == 1
%     posu = strfind(name,'_');
%     poss = posu(find(posu < pos,1,'last'));
%     nFramesStr = name((poss+1):(pos+1));
%     pos = strfind(lower(nFramesStr),'fr');
%     props.number_of_frames = str2double(nFramesStr(1:(pos(1)-1)));
% else
% %     display(name);
% %     xx = input('Enter number of frames = ');
%     props.number_of_frames = -1;
% end

posu = strfind(name,'_');
ps = 1;
pe = posu(1)-1;
props{1} = name(ps:pe);
for ii = 1:length(posu)
    ps = posu(ii)+1;
    if ii < length(posu)
        pe = posu(ii+1)-1;
    else
        pe = length(name);
    end
    props{ii+1} = name(ps:pe);
end

