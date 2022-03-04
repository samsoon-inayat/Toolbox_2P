function rasters = get_rasters_data(ei,selContexts,rasterNames)
if ~exist('ei','var')
    ei = evalin('base','d15_2');
    selContexts = [1 2 3 3 4 4 5 5 6 7];
    rasterNames = {'light22T','air55T','air77T','airD','air77T','airD','air77T','airD','light22T','air55T'};
end

rasters = cell(length(ei),length(selContexts));
for ii = 1:length(ei)
    if sum(selContexts)==0
        rasters(ii,:) = get_data_motion(ei{ii},selContexts,rasterNames);
    else
        rasters(ii,:) = get_data(ei{ii},selContexts,rasterNames);
    end
end

function rasters = get_data_motion(ei,selContexts,rasterNames)
nplanes = length(ei.plane);
for ii = 1:length(selContexts)
    pp = 1;
    if strcmp(rasterNames{ii},'motionOnsets')
        tempR = ei.plane{pp}.motionOnset_rasters;
    else
        tempR = ei.plane{pp}.motionOffset_rasters;
    end
    iscell1 = ei.plane{pp}.tP.iscell(:,1);
    if nplanes == 1
        rasters{ii,1} = tempR;
        rasters{ii,1}.iscell = iscell1;
    else
        pp = 2;
        if strcmp(rasterNames{ii},'motionOnsets')
            tempR1 = ei.plane{pp}.motionOnset_rasters;
        else
            tempR1 = ei.plane{pp}.motionOffset_rasters;
        end
        rasters{ii,1} = combine_planes_data(tempR,tempR1);
        iscell2 = ei.plane{pp}.tP.iscell(:,1);
        rasters{ii,1}.iscell = cat(1,iscell1,iscell2);
        rasters{ii,1}.cns = [[1:length(rasters{ii,1}.iscell)]' [ones(length(iscell1),1);2*ones(length(iscell2),1)] [[1:length(iscell1)]';[1:length(iscell2)]']];
    end
    rasters{ii,1}.context_info = sprintf('%d-%s',selContexts(ii),rasterNames{ii});
    rasters{ii,1}.marker_name = sprintf('%s',rasterNames{ii});
    rasters{ii,1}.thorexp = ei.thorExp;
    rasters{ii,1}.beltLength = get_belt_length(ei);
    spR = rasters{ii,1}.sp_rasters1; mspR = squeeze(nanmean(spR,1));
    [valsT,valsTi] = (max(spR,[],2));
    [vals,valsi] = max(mspR);
    rasters{ii,1}.peak_location = valsi * rasters{ii,1}.bin_width;
    rasters{ii,1}.peak_location_trials = squeeze(valsTi) * rasters{ii,1}.bin_width;
end

function rasters = get_data(ei,selContexts,rasterNames)

nplanes = length(ei.plane);
for ii = 1:length(selContexts)
    if selContexts(ii) == 0
        rasters(ii,1) = get_data_motion(ei,selContexts(ii),rasterNames(ii));
        continue;
    end
    contextNumber = selContexts(ii);
    pp = 1;
    thisContext = ei.plane{pp}.contexts(selContexts(ii));
    disp(sprintf('%s - plane-%d',thisContext.name,pp));
    cmdTxt = sprintf('tempR = thisContext.rasters.%s;',rasterNames{ii});
    eval(cmdTxt);
    if isempty(tempR)
        thisStimMarker = rasterNames{ii}(1:(end-1));
        cmdTxt = sprintf('markersOn = thisContext.markers.%s_onsets;',thisStimMarker);
        eval(cmdTxt);
        cmdTxt = sprintf('markersOff = thisContext.markers.%s_offsets;',thisStimMarker);
        eval(cmdTxt);
        binwidths = evalin('base','binwidths');
        if strcmp(rasterNames{ii}(end),'T')
            thisRasterType = 'time';
        end
        if strcmp(rasterNames{ii}(end),'D')
            thisRasterType = 'dist';
        end
        tempR = make_rasters_quick(ei,pp,markersOn,markersOff,thisRasterType,binwidths);
        trials = 1:size(tempR.sp_rasters,1);
        if strcmp(rasterNames{ii}(end),'T')
            thispFolder = ei.plane{pp}.folder;
            tempR = findRasterProperties_1(thispFolder,contextNumber,thisStimMarker,tempR,thisRasterType,trials,[0 0 0]);
        end
        if strcmp(rasterNames{ii}(end),'D')
             if isfield(ei.plane{pp},'folderD')
                thispFolderD = ei.plane{pp}.folderD;
                if ~exist(thispFolderD)
                    mkdir(thispFolderD);
                end
            else
                thispFolderD = ei.plane{pp}.folder;
            end
            tempR = findRasterProperties_1(thispFolderD,contextNumber,thisStimMarker,tempR,thisRasterType,trials,[0 0 0]);
        end
        
    end
    iscell1 = ei.plane{pp}.tP.iscell(:,1);
    if nplanes == 1
        rasters{ii,1} = tempR;
        rasters{ii,1}.iscell = iscell1;
    else
        pp = 2;
        thisContext = ei.plane{pp}.contexts(selContexts(ii));
        disp(sprintf('%s - plane-%d',thisContext.name,pp));
        cmdTxt = sprintf('tempR1 = thisContext.rasters.%s;',rasterNames{ii});
        eval(cmdTxt);
        if isempty(tempR1)
%             tempR1 = make_rasters(ei,2,markersOn,markersOff,thisRasterType,binwidths);
            cmdTxt = sprintf('markersOn = thisContext.markers.%s_onsets;',thisStimMarker);
            eval(cmdTxt);
            cmdTxt = sprintf('markersOff = thisContext.markers.%s_offsets;',thisStimMarker);
            eval(cmdTxt);
            binwidths = evalin('base','binwidths');
            if strcmp(rasterNames{ii}(end),'T')
                thisRasterType = 'time';
            end
            if strcmp(rasterNames{ii}(end),'D')
                thisRasterType = 'dist';
            end
            tempR1 = make_rasters_quick(ei,pp,markersOn,markersOff,thisRasterType,binwidths);
            trials = 1:size(tempR1.sp_rasters,1);
            if strcmp(rasterNames{ii}(end),'T')
                thispFolder = ei.plane{pp}.folder;
                tempR1 = findRasterProperties_1(thispFolder,contextNumber,thisStimMarker,tempR1,thisRasterType,trials,[0 0 0]);
            end
            if strcmp(rasterNames{ii}(end),'D')
                 if isfield(ei.plane{pp},'folderD')
                    thispFolderD = ei.plane{pp}.folderD;
                    if ~exist(thispFolderD)
                        mkdir(thispFolderD);
                    end
                else
                    thispFolderD = ei.plane{pp}.folder;
                end
                tempR1 = findRasterProperties_1(thispFolderD,contextNumber,thisStimMarker,tempR1,thisRasterType,trials,[0 0 0]);
            end
        end
        rasters{ii,1} = combine_planes_data(tempR,tempR1);
        iscell2 = ei.plane{pp}.tP.iscell(:,1);
        rasters{ii,1}.iscell = cat(1,iscell1,iscell2);
        rasters{ii,1}.cns = [[1:length(rasters{ii,1}.iscell)]' [ones(length(iscell1),1);2*ones(length(iscell2),1)] [[1:length(iscell1)]';[1:length(iscell2)]']];
    end
    rasters{ii,1}.context_info = sprintf('%d-%s',selContexts(ii),rasterNames{ii});
    rasters{ii,1}.marker_name = sprintf('%s',rasterNames{ii});
    rasters{ii,1}.thorexp = ei.thorExp;
    rasters{ii,1}.beltLength = get_belt_length(ei);
    spR = rasters{ii,1}.sp_rasters1; mspR = squeeze(nanmean(spR,1));
    [valsT,valsTi] = (max(spR,[],2));
    [vals,valsi] = max(mspR);
    rasters{ii,1}.peak_location = valsi * rasters{ii,1}.bin_width;
    rasters{ii,1}.peak_location_trials = squeeze(valsTi) * rasters{ii,1}.bin_width;
    n = 0;
end

function tempRC = combine_planes_data(tempR,tempR1)
struct_fields_to_combine = {'fromFrames','info_metrics','gauss_fit_on_mean','fractal_dim'};
fields = fieldnames(tempR);
cellmatflag = 1;
for ff = 1:length(fields)
    thisFieldTxt = fields{ff};
    cmdTxt = sprintf('thisField = tempR.%s;',fields{ff});
    eval(cmdTxt);
    cmdTxt = sprintf('thisField1 = tempR1.%s;',fields{ff});
    eval(cmdTxt);
    if isstruct(thisField)
        if strcmp(fields{ff},'wholeData')
            thisFieldC = combine_whole_data(thisField,thisField1);
        end
        if strcmp(fields{ff},'fromFrames')
            thisFieldC = combine_planes_data(thisField,thisField1);
        end
        if strcmp(fields{ff},'info_metrics')
            thisFieldC = combine_info_metrics(thisField,thisField1);
        end
        if strcmp(thisFieldTxt,'gauss_fit_on_mean')
            thisFieldC = combine_gauss_fit(thisField,thisField1);
        end
        if strcmp(thisFieldTxt,'fractal_dim')
            thisFieldC = combine_fractal_dim(thisField,thisField1);
        end
        n = 0;
    else
        if size(thisField,3) == 1
            if strcmp(fields{ff},'activity_speed_corr')
                thisFieldC = [thisField;thisField1];
            else
                thisFieldC = thisField;
            end
        else
            try
                thisFieldC = cat(3,thisField,thisField1);
            catch
                szf1 = size(thisField,2); szf2 = size(thisField1,2);
                df = szf1 - szf2;
                if df == 0
                    error;
                end
                if df > 0
                    thisFieldC = cat(3,thisField(:,1:szf2,:),thisField1);
                end
                if df < 0
                    thisFieldC = cat(3,thisField,thisField1(:,1:szf1,:));
                end
            end
        end
        if cellmatflag
            cell_numbers = [[1:size(thisFieldC,3)]' [ones(size(thisField,3),1);(2*ones(size(thisField1,3),1))] [[1:size(thisField,3)]';[1:size(thisField1,3)]']];
            cellmatflag = 0;
        end
    end
    cmdTxt = sprintf('tempRC.%s = thisFieldC;',fields{ff});
    eval(cmdTxt);
end
tempRC.cell_numbers = cell_numbers;

function ifC = combine_info_metrics(if1,if2)
% fields = fieldnames(if1);
% for ff = 1:length(fields)
%     cmdTxt = sprintf('ifC.%s = [];',fields{ff}); eval(cmdTxt);
% end
ifC.ShannonMI = cat(2,if1.ShannonMI,if2.ShannonMI);
ifC.ShannonMI_p = cat(2,if1.ShannonMI_p,if2.ShannonMI_p);
ifC.ShannonMI_Zsh = cat(2,if1.ShannonMI_Zsh,if2.ShannonMI_Zsh);
ifC.ShannonMI_shuffle = cat(1,if1.ShannonMI_shuffle,if2.ShannonMI_shuffle);


function thisFieldC = combine_gauss_fit(thisField,thisField1)
thisFieldC.coefficients_Rs_mean = cat(1,thisField.coefficients_Rs_mean,thisField1.coefficients_Rs_mean);
thisFieldC.worked = cat(1,thisField.worked,thisField1.worked);
thisFieldC.gauss1Formula = thisField.gauss1Formula;
if isfield(thisField,'coefficients_Rs_trials')
    thisFieldC.coefficients_Rs_trials = cat(3,thisField.coefficients_Rs_trials,thisField1.coefficients_Rs_trials);
end

function thisFieldC = combine_fractal_dim(thisField,thisField1)
thisFieldC.HaFD = cat(2,thisField.HaFD,thisField1.HaFD);
thisFieldC.HiFD = cat(2,thisField.HiFD,thisField1.HiFD);


function thisFieldC = combine_whole_data(thisField,thisField1)
fields = fieldnames(thisField);
for ii = 1:length(fields)
    fieldTxt = fields{ii};
    if strcmp(fieldTxt,'sp_rasters')
        thisFieldC.sp_rasters = cat(2,thisField.sp_rasters,thisField1.sp_rasters);
    elseif strcmp(fieldTxt,'cell_history')
        thisFieldC.cell_history = cat(1,thisField.cell_history,thisField1.cell_history);
    else
        cmdTxt = sprintf('thisFieldC.%s = thisField.%s;',fields{ii},fields{ii});
        eval(cmdTxt);
    end
end

