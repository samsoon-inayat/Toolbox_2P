function oei = combineRecordings(ei)
cellsInds = [];
cellsIndsRecs = [];
for ii = 1:length(ei)
    cellsN(ii) = size(ei{ii}.signals,1);
    cellsInds = [cellsInds 1:cellsN(ii)];
    cellsIndsRecs = [cellsIndsRecs ii*ones(1,cellsN(ii))];
end

if length(ei) == 1
    oei = ei{1};
    oei.cellsInds = [1:cellsN(1)];
    oei.cellsIndsRecs = ones(1,cellsN(1));
    oei.cellsN = cellsN;
    return;
end

n = 0;



for ii = 1:length(ei)
    tei = ei{ii};
    contexts{ii} = tei.contexts;
end

rTypes = {'D','T'};

for ii = 1:length(contexts{1}) 
    thisContext1 = contexts{1}(ii);
    for mm = 2:length(ei)
        thisContext2 = contexts{mm}(ii);
        stimMarkers = thisContext1.stimMarkers;
        trials = thisContext1.trials;
        if isempty(trials)
            continue;
        end
        for jj = 1:length(stimMarkers)
            if isequal([ii mm jj],[3 2 1])
                n = 0;
            end
            thisStimMarker = stimMarkers{jj};
            if strcmp(thisStimMarker,'air') || strcmp(thisStimMarker,'belt') || strcmp(thisStimMarker,'light')% || strcmp(thisStimMarker,'motionI') || strcmp(thisStimMarker,'motionOnsetsOffsets')
                for kk = 1:length(rTypes)
                    try
                        cmdTxt = sprintf('temp1 = thisContext1.rasters.%s%s;',stimMarkers{jj},rTypes{kk});
                        eval(cmdTxt);
                    catch
                        continue;
                    end
                    cmdTxt = sprintf('temp2 = thisContext2.rasters.%s%s;',stimMarkers{jj},rTypes{kk});
                    eval(cmdTxt);
                    cTemp = temp1;
                    [rows1,cols1] = size(temp1.rasters(:,:,1));
                    [rows2,cols2] = size(temp2.rasters(:,:,1));
                    rows = min(rows1,rows2);
                    cols = min(cols1,cols2);
                    cTemp.rasters = cat(3,temp1.rasters(1:rows,1:cols,:),temp2.rasters(1:rows,1:cols,:));
                    cTemp.SI = cat(2,temp1.SI,temp2.SI);
                    if isfield(cTemp,'gauss_fit_on_mean')
                        cTemp.gauss_fit_on_mean.allgof = cat(2,temp1.gauss_fit_on_mean.allgof,temp2.gauss_fit_on_mean.allgof);
                        cTemp.gauss_fit_on_mean.alloutput = cat(2,temp1.gauss_fit_on_mean.alloutput,temp2.gauss_fit_on_mean.alloutput);
                        cTemp.gauss_fit_on_mean.coeffs = cat(1,temp1.gauss_fit_on_mean.coeffs,temp2.gauss_fit_on_mean.coeffs);
                        cTemp.gauss_fit_on_mean.worked = cat(2,temp1.gauss_fit_on_mean.worked,temp2.gauss_fit_on_mean.worked);
                    end
                    cTemp.belt_length = ei{1}.b.belt_length;
                    cmdTxt = sprintf('thisContext1.rasters.%s%s = cTemp;',stimMarkers{jj},rTypes{kk});
                    eval(cmdTxt);
                end
            end
            if strcmp(thisStimMarker,'motionOnsets22') || strcmp(thisStimMarker,'motionOffsets22') || strcmp(thisStimMarker,'airOnsets22') ...
                    || strcmp(thisStimMarker,'airOffsets22') || strcmp(thisStimMarker,'airI') || strcmp(thisStimMarker,'motionOffsetAirOnset')
                for kk = 2
                    try
                        cmdTxt = sprintf('temp1 = thisContext1.rasters.%s%s;',stimMarkers{jj},rTypes{kk});
                        eval(cmdTxt);
                    catch
                        continue;
                    end
                    cmdTxt = sprintf('temp2 = thisContext2.rasters.%s%s;',stimMarkers{jj},rTypes{kk});
                    eval(cmdTxt);
                    cTemp = temp1;

                    [rows1,cols1] = size(temp1.rasters(:,:,1));
                    [rows2,cols2] = size(temp2.rasters(:,:,1));
                    rows = min(rows1,rows2);
                    cols = min(cols1,cols2);
                    cTemp.rasters = cat(3,temp1.rasters(1:rows,1:cols,:),temp2.rasters(1:rows,1:cols,:));
                    cTemp.SI = cat(2,temp1.SI,temp2.SI);
                    if isfield(cTemp,'gauss_fit_on_mean')
                        cTemp.gauss_fit_on_mean.allgof = cat(2,temp1.gauss_fit_on_mean.allgof,temp2.gauss_fit_on_mean.allgof);
                        cTemp.gauss_fit_on_mean.alloutput = cat(2,temp1.gauss_fit_on_mean.alloutput,temp2.gauss_fit_on_mean.alloutput);
                        cTemp.gauss_fit_on_mean.coeffs = cat(1,temp1.gauss_fit_on_mean.coeffs,temp2.gauss_fit_on_mean.coeffs);
                        cTemp.gauss_fit_on_mean.worked = cat(2,temp1.gauss_fit_on_mean.worked,temp2.gauss_fit_on_mean.worked);
                    end
                    cmdTxt = sprintf('thisContext1.rasters.%s%s = cTemp;',stimMarkers{jj},rTypes{kk});
                    eval(cmdTxt);
                end
            end
        end
    end
    contexts{1}(ii) = thisContext1;
end
% oei.eis = ei;
oei.contexts = contexts{1};
oei.cellsInds = cellsInds;
oei.cellsIndsRecs = cellsIndsRecs;
oei.cellsN = cellsN;
n = 0;