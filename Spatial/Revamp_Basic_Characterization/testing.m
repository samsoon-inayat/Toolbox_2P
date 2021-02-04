function motion_responsive_cells(ei)

if ~exist('ei','var')
    ei = evalin('base','ei1(1)');
end

cai_sampling_rate = ei{1}.thorExp.frameRate;
effective_sampling_rate = cai_sampling_rate;
% effective_sampling_rate = 1/0.1;

pp = 1;

condN = 10;
rasters_ao = ei{1}.plane{pp}.contexts(condN).rasters.air33T.sp_rasters1;%    fromFrames.sp_rasters;

condN = 9;
rasters = ei{1}.plane{pp}.contexts(condN).rasters.light11T.sp_rasters1;%    fromFrames.sp_rasters;
rasters = ei{1}.plane{pp}.contexts(condN).rasters.light11T.fromFrames.sp_rasters;%    fromFrames.sp_rasters;

condN = 1;
rastersD1 = ei{1}.plane{pp}.contexts(condN).rasters.airD.sp_rasters1;%    fromFrames.sp_rasters;
zMIs1 = ei{1}.plane{pp}.contexts(condN).rasters.airD.MIs;

condN = 3;
rastersD2 = ei{1}.plane{pp}.contexts(condN).rasters.airD.sp_rasters1;%    fromFrames.sp_rasters;
zMIs2 = ei{1}.plane{pp}.contexts(condN).rasters.airD.MIs;

condN = 6;
rastersD3 = ei{1}.plane{pp}.contexts(condN).rasters.airD.sp_rasters1;%    fromFrames.sp_rasters;
zMIs3 = ei{1}.plane{pp}.contexts(condN).rasters.airD.MIs;

timeBefore = 1;
number_of_columns = size(rasters,2);
column_index = round(timeBefore * effective_sampling_rate);
cis = [];
cis(1,:) = [1 column_index (number_of_columns-column_index+1)];
column_index = cis - 1; column_index(1) = []; column_index(3) = number_of_columns;
cis(2,:) = column_index;
resp = findResponsiveRasters(rasters,cis);
plotRasters_multi({rastersD1,rastersD2,rastersD3,rasters},find(resp.p < 0.05),[])
plotRasters_multi({rastersD1,rastersD2,rastersD3,rasters},find(zMIs3'>0.2),[])
return;
% props = raster_properties(rasters);
% clus = getClusters(props);
n = 0;

n_samples = size(rasters,2);
totalTime = n_samples/effective_sampling_rate;
mdata.xs = linspace(0,totalTime,n_samples);
% % mdata.xs = round(mdata.xs - mdata.xs(column_index),1);
mdata.cis = cis(1,2:3);
plotRasters_simple(rasters,find(resp.ps(:,1)<0.05),mdata)
plotRasters_simple(rasters,find(resp.ps(:,2)<0.05),mdata)
plotRasters_multi(rasters,rastersD,find(resp.ps(:,1)<0.05),[])
plotRasters_multi({rastersD1,rastersD2,rastersD3,rasters},find(zMIs1'>5),[])
% plotRasters_simple(rasters,[8 28 19 29 35 47 71 190 197 204 199 186 191 235 282 322 319 344],mdata)

n = 0;