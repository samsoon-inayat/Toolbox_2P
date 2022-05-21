function ei = make_and_load_rasters(ei,binwidths,owr,contextsIncExc)

if ~exist('owr','var')
    owr = [0 0 0];
end
if ~exist('contextsIncExc','var')
    contextsIE_val = 0;
    stimMarkerIE_val = 0;
else
    contextsIE_val = contextsIncExc{1};         contextsIE = contextsIncExc{2};
    if length(contextsIncExc) > 2
        stimMarkerIE_val = contextsIncExc{3};         stimMarkerIE = contextsIncExc{4};
    else
        stimMarkerIE_val = 0;
    end
end

allContexts = contextDefinitions;
for aa = 1:length(ei)
    ptei = ei{aa};
    tei = ptei;
    for pp = 1:length(ei{aa}.plane)
        thispFolder = tei.plane{pp}.folder;
        if isfield(tei.plane{pp},'folderD')
            thispFolderD = tei.plane{pp}.folderD;
            if ~exist(thispFolderD)
                mkdir(thispFolderD);
            end
        else
            thispFolderD = tei.plane{pp}.folder;
        end
        disp(thispFolder);
        
        tei.folders.thispFolder = thispFolder;
%         tei.areCells = find(tei.plane{pp}.tP.areCells);
%         tei.deconv = tei.plane{pp}.tP.deconv;
%         tei.ops1 = tei.plane{pp}.tP.ops;
%         tei.tP = tei.plane{pp}.tP;
        tplane = tei.plane{pp};
        contexts = getContexts({ptei},'define_cont_11.m');
        for contextNumber = 1:length(contexts)
            thisContext = contexts(contextNumber);
            if contextsIE_val == -1
                if found(thisContext.name,contextsIE)
                    continue;
                end
            end
            if contextsIE_val == 1
                if ~found(thisContext.name,contextsIE)
                    continue;
                end
            end
            stimMarkers = thisContext.stimMarkers;
            trials = thisContext.trials;
            if isempty(trials)
                continue;
            end
            if length(stimMarkers) == 0
                continue;
            end
            for jj = 1:length(stimMarkers)
                thisStimMarker = stimMarkers{jj};
                if stimMarkerIE_val == -1
                    if exact_match(thisStimMarker,stimMarkerIE)
                        continue;
                    end
                end
                if stimMarkerIE_val == 1
                    if ~exact_match(thisStimMarker,stimMarkerIE)
                        continue;
                    end
                end
             
                cmdTxt = sprintf('markersOn = contexts(contextNumber).markers.%s_onsets;',thisStimMarker);
                eval(cmdTxt);
                cmdTxt = sprintf('markersOff = contexts(contextNumber).markers.%s_offsets;',thisStimMarker);
                eval(cmdTxt);
                if isempty(markersOn) & isempty(markersOff)
                    continue;
                end
                typesOfRasters = getTypesOfRasters(allContexts,thisStimMarker);
                for kk = 1:length(typesOfRasters)
                    thisRasterType = typesOfRasters{kk};
                    disp(sprintf('%s--- %s -- %s',thisContext.name,stimMarkers{jj},thisRasterType));
                    rasters = make_rasters(tei,pp,markersOn,markersOff,thisRasterType,binwidths);
                    trials = 1:size(rasters.sp_rasters,1);
                    if length(unique(rasters.lastBin)) > 1
                        n= 0;
                    end
                    if strcmp(thisRasterType,'dist')
                        rasters = findRasterProperties_1(thispFolderD,contextNumber,thisStimMarker,rasters,thisRasterType,trials,owr);
                        if double(rasters.bin_width) == double(binwidths(2))
                        else
                            error;
                        end
                        cmdTxt = sprintf('contexts(contextNumber).rasters.%sD = rasters;',thisStimMarker); eval(cmdTxt);
                    end
                   if strcmp(thisRasterType,'time')
%                        if strcmp(thisStimMarker,'airI')
%                            owr(1) = 1;
%                        else
%                            owr(1) = 0;
%                        end
                        rasters = findRasterProperties_1(thispFolder,contextNumber,thisStimMarker,rasters,thisRasterType,trials,owr);
                        if double(rasters.bin_width) == double(binwidths(1))
                        else
                            error;
                        end
                        cmdTxt = sprintf('contexts(contextNumber).rasters.%sT = rasters;',thisStimMarker); eval(cmdTxt);
                    end 
                end
            end
        end
        ptei.plane{pp}.contexts = contexts;
%         ptei.b.belt_length = temp.belt_length;
    end
    ei{aa} = ptei;
end



function fo = found(name,contextsIE)
fo = 0;
for cniii = 1:length(contextsIE)
    cniiind = strfind(name,contextsIE{cniii});
    if ~isempty(cniiind)
        fo = 1;
        break;
    end
end

function fo = exact_match(name,contextsIE)
fo = 0;
if sum(strcmp(contextsIE,name)) > 0
    fo = 1;
end