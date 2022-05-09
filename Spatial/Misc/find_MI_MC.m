function ei = find_MI_MC(ei,owr,contextsIncExc)

if ~exist('ei','var')
    ei = evalin('base','ei');
end

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
        else
            thispFolderD = tei.plane{pp}.folder;
        end
        disp(thispFolder);
        
        tei.folders.thispFolder = thispFolder;
        tplane = tei.plane{pp};
        contexts = tplane.contexts; contextsMC = tplane.contextsMC;
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
            thisContextMC = contextsMC(contextNumber);
            stimMarkers = thisContext.stimMarkers;
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
                typesOfRasters = getTypesOfRasters(allContexts,thisStimMarker);
                for kk = 1:length(typesOfRasters)
                    thisRasterType = typesOfRasters{kk};
                    disp(sprintf('%s--- %s -- %s',thisContext.name,stimMarkers{jj},thisRasterType));
                    if strcmp(thisRasterType,'dist')
                        cmdTxt = sprintf('rasters = contexts(contextNumber).rasters.%sD;',thisStimMarker); eval(cmdTxt);
                        cmdTxt = sprintf('rastersMC = contextsMC(contextNumber).rasters.%sD;',thisStimMarker); eval(cmdTxt);
                    end
                   if strcmp(thisRasterType,'time')
                        cmdTxt = sprintf('rasters = contexts(contextNumber).rasters.%sT;',thisStimMarker); eval(cmdTxt);
                        cmdTxt = sprintf('rastersMC = contextsMC(contextNumber).rasters.%sT;',thisStimMarker); eval(cmdTxt);
                   end 
                   trials = 1:size(rasters.sp_rasters,1);
                   rasters = findRasterProperties_1_MC(thispFolder,contextNumber,thisStimMarker,{rasters,rastersMC},thisRasterType,trials,owr);
                   if strcmp(thisRasterType,'dist')
                        cmdTxt = sprintf('contexts(contextNumber).rasters.%sD = rasters;',thisStimMarker); eval(cmdTxt);
                    end
                   if strcmp(thisRasterType,'time')
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
