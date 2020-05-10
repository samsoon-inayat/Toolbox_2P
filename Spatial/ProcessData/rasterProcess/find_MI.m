function out =  find_MI(thispFolder,contextNumber,stimMarker,rasters,thisRasterType,trials,ow)

fileName = makeName(sprintf('info_metrics_%d_%s_%s.mat',contextNumber,stimMarker,thisRasterType),thispFolder);
if exist(fileName,'file') && ow == 0
    out = load(fileName);
    return;
end
if isfield(rasters,'sp_rasters_nan_corrected')
    Rs = rasters.sp_rasters_nan_corrected(trials,:,:);
    Dur = rasters.duration_nan_corrected(trials,:,:);
else
    out = [];
    return;
end

nbins = 4;
nShuffle = 500;
hWaitBar = waitbar(0,sprintf('Finding MI scores ... plz wait'));
numCells = size(Rs,3);

[output,~] = info_metrics_S_onlyMI(Rs(:,:,1),[],nbins,Dur,nShuffle);
fields = fieldnames(output);

O(numCells,1) = struct;

for ff = 1:length(fields)
    thisField = fields{ff};
    cmdTxt = sprintf('O(1).%s = [];',thisField);
    eval(cmdTxt);
end

tic
parfor ii = 1:numCells
    rng(3,'twister');
%     waitbar(ii/numCells,hWaitBar,sprintf('MI - Processing cell %d/%d',ii,numCells));
    [O(ii),~] = info_metrics_S_onlyMI(Rs(:,:,ii),[],nbins,Dur,nShuffle);
end
close(hWaitBar);
toc

for ii = 1:numCells
    output = O(ii);
    tempVar = [];
    for ff = 1:length(fields)
        thisField = fields{ff};
        cmdTxt = sprintf('tempVar = output.%s;',thisField);
        eval(cmdTxt);
        if length(tempVar) == 1
            cmdTxt = sprintf('out.%s(ii) = tempVar;',thisField);
            eval(cmdTxt);
        else
            cmdTxt = sprintf('out.%s(ii,:) = tempVar;',thisField);
            eval(cmdTxt);
        end
    end
end

save(fileName,'-struct','out','-v7.3');


