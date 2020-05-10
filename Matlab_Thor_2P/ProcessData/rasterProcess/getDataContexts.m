function [ddata,tdata] = getDataContexts(ei,selContexts,markerType)

for ii = 1:length(selContexts)
%     ii
    thisContext = ei.contexts(selContexts(ii));
    stimMarkers = thisContext.stimMarkers;
    trials = thisContext.trials;
%     for jj = 1:length(stimMarkers)
        thisStimMarker = markerType;%stimMarkers{jj};
        if strcmp(thisStimMarker,'air') || strcmp(thisStimMarker,'belt') || strcmp(thisStimMarker,'motionI') || strcmp(thisStimMarker,'motionOnsetsOffsets')
            cmdTxt = sprintf('ddata{ii} = thisContext.rasters.%sD; ddata{ii}.name = thisContext.name; ddata{ii}.stimMarker = thisStimMarker;',thisStimMarker);
            try
                eval(cmdTxt);
%                 [pcw,pcc] = getPlaceCellProps(ddata{ii},1:size(ddata{ii}.rasters,3));
%                 ddata{ii}.pcw = pcw;
%                 ddata{ii}.pcc = pcc;
            catch
                ddata{ii} = [];
            end
            cmdTxt = sprintf('tdata{ii} = thisContext.rasters.%sT; tdata{ii}.name = thisContext.name; tdata{ii}.stimMarker = thisStimMarker;',thisStimMarker);
            try
                eval(cmdTxt);
            catch
                tdata{ii} = [];
            end
            
        end
        
        if strcmp(thisStimMarker,'motionOnsets22') || strcmp(thisStimMarker,'motionOffsets22') || strcmp(thisStimMarker,'airOnsets22') ...
                || strcmp(thisStimMarker,'airOffsets22') || strcmp(thisStimMarker,'airI') || strcmp(thisStimMarker,'motionOffsetAirOnset') || strcmp(thisStimMarker,'light')...
                || strcmp(thisStimMarker,'tone') || strcmp(thisStimMarker,'airOnsets27') ...
                || strcmp(thisStimMarker,'airOffsets27') || strcmp(thisStimMarker,'airOnsets11') ...
                || strcmp(thisStimMarker,'airOffsets11') || strcmp(thisStimMarker,'airOnsets010') ...
                || strcmp(thisStimMarker,'airOffsets010') || strcmp(thisStimMarker,'airOnsets01') ...
                || strcmp(thisStimMarker,'airOffsets01') 
            ddata{ii} = [];
            cmdTxt = sprintf('tdata{ii} = thisContext.rasters.%sT; tdata{ii}.name = thisContext.name;',thisStimMarker);
            try
                eval(cmdTxt);
            catch
                tdata{ii} = [];
            end
        end
%     end
end
