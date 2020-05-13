function [folders_data] = getFilesData_EV_wrapper (eDataFolder,stimulus_types,folders_data)
eDataFolder
tfolders = dir(eDataFolder);
folders = tfolders(3:end);
for ii = 1:length(folders)
    if folders(ii).isdir
        thisName = folders(ii).name;
        thisFolder = fullfile(eDataFolder,folders(ii).name);
        temp = dir(thisFolder);
        temp1 = temp(3:end);
        hi1 = findFilesContainingKeyword(temp1,'hi3');
        if isempty(hi1)
            folders_data = getFilesData_EV_wrapper (thisFolder,stimulus_types,folders_data)
            continue;
        else
            iii = length(folders_data)+1;
            folders_data(iii).name = thisName;
            folders_data(iii).root_folder = thisFolder;
            folders_data(iii).stimulus = findStimulusType(thisName,stimulus_types);
            folders_data = getFilesData_EV (thisFolder,stimulus_types,folders_data);
        end
    end
end


function stimulus = findStimulusType (name,stimTypes)

signatures = {'iso','Hz','ms','V','mA','fr'};

posu = strfind(name,'_');
ps = 1;
pe = posu(1)-1;
stimXics{1} = name(ps:pe);
for ii = 1:length(posu)
    ps = posu(ii)+1;
    if ii < length(posu)
        pe = posu(ii+1)-1;
    else
        pe = length(name);
    end
    stimXics{ii+1} = name(ps:pe);
end
stimulus.all_props = stimXics;
sMat = zeros(length(signatures),length(stimXics));
for ii = 1:length(stimXics)
    for jj = 1:length(signatures)
        if ~isempty(strfind(stimXics{ii},signatures{jj}))
           sMat(jj,ii) = 1;
        end
    end
end
[r,c] = find(sMat == 1);
for cc = 1:length(c)
    stimulus.props{cc}.type = signatures{r(cc)};
    stimulus.props{cc}.val = sscanf(stimXics{c(cc)},'%g',5);
end
stimMat = zeros(length(stimTypes),length(stimXics));
for ii = 1:length(stimXics)
    for jj = 1:length(stimTypes)
        if ~isempty(strfind(stimXics{ii},stimTypes{jj}))
           stimMat(jj,ii) = 1;
        end
    end
end
if sum(stimMat(:)) == 1
    [r,c] = find(stimMat == 1);
    stimulus.stimulus_type = stimTypes{r};
end
