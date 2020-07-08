function out =  find_MI_1(thispFolder,contextNumber,stimMarker,rasters,thisRasterType,trials,ow)

fileName = makeName(sprintf('info_metrics_1_%d_%s_%s.mat',contextNumber,stimMarker,thisRasterType),thispFolder);
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
nShuffle = 0;
hWaitBar = waitbar(0,sprintf('Finding MI scores ... plz wait'));
numCells = size(Rs,3);

trial_vals = NaN(size(Rs,1),numCells);
ShannonMI = NaN(1,numCells);
ShannonMI_Zsh = NaN(1,numCells);
tic
parfor ii = 1:numCells
    rng(3,'twister');
%     waitbar(ii/numCells,hWaitBar,sprintf('MI - Processing cell %d/%d',ii,numCells));
    [ShannonMI(1,ii),ShannonMI_Zsh(1,ii),~] = info_metrics_S_onlyMI_1(Rs(:,:,ii),[],nbins,Dur,nShuffle);
    trial_vals(:,ii) = find_MI_trials(Rs(:,:,ii),nbins,Dur,nShuffle);
end
close(hWaitBar);
toc

out.ShannonMI = ShannonMI;
out.ShannonMI_trials = trial_vals;

% for ii = 1:numCells
%     output = O(ii);
%     tempVar = [];
%     for ff = 1:length(fields)
%         thisField = fields{ff};
%         cmdTxt = sprintf('tempVar = output.%s;',thisField);
%         eval(cmdTxt);
%         if length(tempVar) == 1
%             cmdTxt = sprintf('out.%s(ii) = tempVar;',thisField);
%             eval(cmdTxt);
%         else
%             cmdTxt = sprintf('out.%s(ii,:) = tempVar;',thisField);
%             eval(cmdTxt);
%         end
%     end
% end

save(fileName,'-struct','out','-v7.3');

function ShannonMI = find_MI_trials(Rs,nbins,Dur,nShuffle)
% trials(size(Rs,1),1) = struct;
ShannonMI = NaN(size(Rs,1),1);
for ii = 1:size(Rs,1)
    [ShannonMI(ii,1),~,~] = info_metrics_S_onlyMI_1(Rs(ii,:),[],nbins,Dur(ii,:),nShuffle);
end
% 
% trial_vals = zeros(length(trials),2);
% for ii = 1:length(trials)
%     trial_vals(ii,1) = trials(ii).ShannonMI;
%     trial_vals(ii,2) = trials(ii).ShannonMI_Zsh;
% end

