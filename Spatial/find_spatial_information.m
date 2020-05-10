function out =  find_spatial_information(ei,rasters,ow)

% out.SIm = find_SI(rasters);
% return;

if ~exist('ow','var')
    ow = 0;
end

if iscell(ei)
    for ii = 1:length(ei)
        N(ii) = length(ei{ii}.areCells);
        if ii == 1
            s = 1;
            e = N(ii);
        else
            s = N(ii-1)+1;
            e = N(ii-1)+N(ii);
        end
        inds = s:e;
        rastersT = rasters;
        rastersT.cells = rastersT.cells(inds);
        out = find_spatial_information(ei{ii},rastersT,ow);
        if ii == 1
            SIs = out.SI;
%             MIs = out.MI;
        else
            SIs = [SIs out.SI];
%             MIs = [MIs out.MI];
        end
    end
    out.SI = SIs;
%     out.MI = MIs;
    return;
end

onsets = rasters.onsets;
offsets = rasters.offsets;
fileName = makeName(sprintf('spatial_information_%d_%d_%d_%d.mat',onsets(1),onsets(end),offsets(1),offsets(end)),ei.folders.thispFolder);
if exist(fileName,'file') && ow == 0
    out = load(fileName);
    return;
end

cellsR = rasters.cells;

for ii = 1:length(cellsR)
    spSignal = ei.deconv.spSigAll{ii};
    caSignal = ei.deconv.caSigAll{ii};
    out.SIm(ii) = spatial_information(ei.b.ts(ei.b.frames_f(1:length(spSignal))),spSignal,cellsR(ii).raster_spikeCount,cellsR(ii).raster_meanFiringRate,rasters.duration,'Dun');
    out.MI(ii) = spatial_information(ei.b.ts(ei.b.frames_f(1:length(spSignal))),spSignal,cellsR(ii).raster_spikeCount,cellsR(ii).raster_meanFiringRate,rasters.duration,'MI_formula');
%     out.SIm(ii) = spatial_informationCa(ei.b.ts(ei.b.frames_f),caSignal,cellsR(ii).raster_spikeCount,cellsR(ii).raster_df,rasters.duration,'Dun');
    lambda(ii,:) = nanmean(cellsR(ii).raster_meanFiringRate./rasters.duration);
    out.M(ii) = max(lambda(ii,:));
end
lambda = lambda';
Pi = mean(rasters.duration);
Pi = repmat(Pi,size(lambda,2),1);


lamb=lambda;
m_lamb=mean(lamb);
Pi=Pi./sum(Pi,2);
Pi=Pi';

SI_series=Pi.*lamb./m_lamb.*log2(lamb./m_lamb);
out.SI=nansum(SI_series);

save(fileName,'-struct','out','-v7.3');

