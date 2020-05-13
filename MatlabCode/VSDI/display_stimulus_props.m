function get_stimulus_props (ags,pags,stim_folders_info,indices,prop_name)
try
if strcmp(prop_name,'All')
    for ii = 1:length(indices)
        ss = stim_folders_info.list(indices(ii),1);
        gg = stim_folders_info.list(indices(ii),2);
        aa = stim_folders_info.list(indices(ii),3);
        rr = stim_folders_info.list(indices(ii),4);
        display(ags(gg).animals{aa}.eRecordings{rr}.stimulus.all_props)
    end
    return;
end

if strcmp(prop_name,'name')
    for ii = 1:length(indices)
        ss = stim_folders_info.list(indices(ii),1);
        gg = stim_folders_info.list(indices(ii),2);
        aa = stim_folders_info.list(indices(ii),3);
        rr = stim_folders_info.list(indices(ii),4);
        display(ags(gg).animals{aa}.eRecordings{rr}.name)
    end
    return;
end

if strcmp(prop_name,'animal_name')
    for ii = 1:length(indices)
        ss = stim_folders_info.list(indices(ii),1);
        gg = stim_folders_info.list(indices(ii),2);
        aa = stim_folders_info.list(indices(ii),3);
        rr = stim_folders_info.list(indices(ii),4);
        display(ags(gg).animals{aa}.name)
    end
    return;
end

list = [];
for ii = 1:length(indices)
    ss = stim_folders_info.list(indices(ii),1);
    gg = stim_folders_info.list(indices(ii),2);
    aa = stim_folders_info.list(indices(ii),3);
    rr = stim_folders_info.list(indices(ii),4);
    for jj = 1:length(ags(gg).animals{aa}.eRecordings{rr}.stimulus.props)
        propType = ags(gg).animals{aa}.eRecordings{rr}.stimulus.props{jj}.type;
        propVal = ags(gg).animals{aa}.eRecordings{rr}.stimulus.props{jj}.val;
        if strcmp(propType,prop_name)
            list = [list propVal];
            break;
        end
    end
end
list
return;

catch
    display('Unrecognized Property')
end