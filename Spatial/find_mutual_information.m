function out =  find_mutual_information(ei,rasters,ow)

% out.SIm = find_SI(rasters);
% return;

if ~exist('ow','var')
    ow = 0;
end

% if iscell(ei)
%     for ii = 1:length(ei)
%         N(ii) = length(ei{ii}.areCells);
%         if ii == 1
%             s = 1;
%             e = N(ii);
%         else
%             s = N(ii-1)+1;
%             e = N(ii-1)+N(ii);
%         end
%         inds = s:e;
%         rastersT = rasters;
%         rastersT.cells = rastersT.cells(inds);
%         out = find_mutual_information(ei{ii},rastersT,ow);
%         if ii == 1
%             SIs = out.SI;
%             MIs = out.MI;
% %             sMIs = out.sMI;
%             zMIs = out.zMI;
%         else
%             SIs = [SIs out.SI];
%             MIs = [MIs out.MI];
% %             sMIs = [sMIs;out.sMI];
%             zMIs = [zMIs out.zMI];
%         end
%     end
%     out.SI = SIs;
%     out.MI = MIs;
% %     out.sMI = sMIs;
%     out.zMI = zMIs;
%     return;
% end

onsets = rasters.onsets;
offsets = rasters.offsets;
fileName = makeName(sprintf('mutual_information_%d_%d_%d_%d.mat',onsets(1),onsets(end),offsets(1),offsets(end)),ei.folders.thispFolder);
if exist(fileName,'file') && ow == 0
    out = load(fileName);
    return;
end
cellsR = rasters.cells;

for ii = 1:length(cellsR)
    if isfield(cellsR(ii),'raster_meanFiringRate')
        allRs(:,:,ii) = cellsR(ii).raster_meanFiringRate./rasters.duration;
    else
        allRs(:,:,ii) = cellsR(ii).sigRaster;
    end
end
allRs = reshape(allRs,1,numel(allRs));
vals = quantile(allRs(~isnan(allRs)),[0.25 0.5 0.75 1]);

maxVal = quantile(allRs(~isnan(allRs)),0.99);

nbins = 4;
FR_threshold = 0.3;
FR_increment = 0.5;
FR_bin_s = linspace(FR_threshold,maxVal,nbins);
FR_bin_e = linspace(FR_increment,maxVal,nbins);

% FR_bin_s = FR_threshold:FR_increment:maxVal;
% FR_bin_e = FR_increment:FR_increment:maxVal;
% FR_bin_s = FR_bin_s(1:nbins);
% FR_bin_e = FR_bin_e(1:nbins);
% Pi = mean(rasters.duration);
% Pi=Pi./sum(Pi,2);
hWaitBar = waitbar(0,sprintf('Finding MI scores ... plz wait'));
MI = [];
zMI = [];
for ii = 1:length(cellsR)
    if ii == 21
        n = 0;
    end
    waitbar(ii/length(cellsR),hWaitBar,sprintf('MI - Processing cell %d/%d',ii,length(cellsR)));
%     thisRaster = cellsR(ii).raster_meanFiringRate./rasters.duration;
%     [out.MI(ii),out.sMI(ii,:),out.zMI(ii)] = findMI(cellsR(ii).raster_meanFiringRate,rasters.duration,nbins,FR_bin_s,FR_bin_e);
    rng(3,'twister');
    if isfield(cellsR(ii),'raster_meanFiringRate')
%         [output,FRbin] =
%         info_metrics_S(cellsR(ii).raster_meanFiringRate./rasters.duration,[],nbins,rasters.duration,100);
%         %changed January 26, 2020 ... not dividing by rasters.duration
%         because it is already firing rate
        [output,FRbin] = info_metrics_S(cellsR(ii).raster_meanFiringRate,[],nbins,rasters.duration,100);
%         [output,FRbin] = info_metrics_S(cellsR(ii).raster_meanFiringRate(1:10,:),[],nbins,rasters.duration(1:10,:),100);
    else
        [output,FRbin] = info_metrics_S(cellsR(ii).sigRaster,[],nbins,rasters.duration,100);
    end
    MI(ii) = output.ShannonMI;
    zMI(ii) = output.ShannonMI_Zsh;
%     [out.MI(ii),~,out.zMI(ii)] = findMI(cellsR(ii).raster_meanFiringRate,rasters.duration,nbins,FR_bin_s,FR_bin_e);
%     
%     figure(1000);clf;plot(nanmean(thisRaster));
%     figure(1000);clf;imagesc(thisRaster);
end
close(hWaitBar);
out.MI = MI;
out.zMI = zMI;

% for ii = 1:length(cellsR)
% %     spSignal = ei.deconv.spSigAll{ii};
% %     caSignal = ei.deconv.caSigAll{ii};
%     if isfield(cellsR(ii),'raster_meanFiringRate')
%         lambda(ii,:) = nanmean(cellsR(ii).raster_meanFiringRate./rasters.duration);
%     else
%         lambda(ii,:) = nanmean(cellsR(ii).sigRaster);
%     end
% %     out.M(ii) = max(lambda(ii,:));
% end
% lambda = lambda';
% Pi = mean(rasters.duration);
% Pi = repmat(Pi,size(lambda,2),1);
% 
% 
% lamb=lambda;
% m_lamb=mean(lamb);
% Pi=Pi./sum(Pi,2);
% Pi=Pi';
% 
% SI_series=Pi.*lamb./m_lamb.*log2(lamb./m_lamb);
% out.SI=nansum(SI_series);

% save(fileName,'-struct','out','-v7.3');

% function [MI sMI zMI] = findMI(thisSig,thisDur,nbins,FR_bin_s,FR_bin_e,varargin)
% thisRaster = (thisSig./thisDur);
% for jj = 1:nbins
%     if jj < nbins
%         temp = thisRaster >= FR_bin_s(jj) & thisRaster < FR_bin_e(jj);
%     else
%         temp = thisRaster >= FR_bin_s(jj);
%     end
%     pBF(jj,:) = sum(temp)/(size(thisRaster,1)*size(thisRaster,2));
% end
% pB = sum(pBF);
% pF = sum(pBF,2);
% for kk = 1:length(pF)
%     for mm = 1:length(pB)
%         sTemp(kk,mm) = pBF(kk,mm)*log2(pBF(kk,mm)/(pF(kk).*pB(mm)));
%     end
% end
% MI = nansum(nansum(sTemp));
% if nargin == 5
%     for aa = 1:100
%         sThisSig = thisSig(randperm(size(thisRaster,1)),randperm(size(thisRaster,2)));
%         sThisDur = thisDur(randperm(size(thisRaster,1)),randperm(size(thisRaster,2)));
% %         sThisRaster = thisRaster(randperm(size(thisRaster,1)),randperm(size(thisRaster,2)));
%         [sMI(aa),~,~] = findMI(sThisSig,sThisDur,nbins,FR_bin_s,FR_bin_e,1);
%     end
%     zMI = normcdf(MI,mean(sMI),std(sMI));
% %     zMI = (MI - mean(sMI))/std(sMI);
% else
%     sMI = NaN;
%     zMI = NaN;
% end