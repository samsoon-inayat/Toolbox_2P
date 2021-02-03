function motion_responsive_cells(ei)

if ~exist('ei','var')
    ei = evalin('base','data(4)');
end

contexts = getContexts(ei);
condition_number = 4;
b = ei{1}.b;
onsets = contexts(condition_number).markers.air_onsets;
offsets = contexts(condition_number).markers.air_offsets;

% onsets = contexts(1).markers.motionI_onsets;
% offsets = contexts(1).markers.motionI_offsets;

% plotMarkers(b,onsets,offsets,100,'motion')
timeBefore = 3; timeAfter = 3;
markersOn = onsets - round(1e6 * timeBefore/b.si);
markersOff = offsets + round(1e6 * timeAfter/b.si);

cai_sampling_rate = ei{1}.thorExp.frameRate;
effective_sampling_rate = cai_sampling_rate;

pp = 1;
ei{1}.tP = ei{1}.plane{pp}.tP;
ei{1}.deconv = ei{1}.plane{pp}.tP.deconv;
raster_data =  getRasters_fixed_bin_width_ctrl(ei{1},pp,markersOn,markersOff,'time');
% raster_data_d =  getRasters_fixed_bin_width_ctrl(ei{1},pp,onsets,offsets,'dist');
rasters = raster_data.fromFrames.sp_rasters;

condN = 1;
rastersD1 = ei{1}.plane{pp}.contexts(condN).rasters.airD.sp_rasters1;%    fromFrames.sp_rasters;
zMIs1 = ei{1}.plane{pp}.contexts(condN).rasters.airD.info_metrics.ShannonMI_Zsh;

condN = 2;
rastersD2 = ei{1}.plane{pp}.contexts(condN).rasters.airD.sp_rasters1;%    fromFrames.sp_rasters;
zMIs2 = ei{1}.plane{pp}.contexts(condN).rasters.airD.info_metrics.ShannonMI_Zsh;

condN = 3;
rastersD3 = ei{1}.plane{pp}.contexts(condN).rasters.airD.sp_rasters1;%    fromFrames.sp_rasters;
zMIs3 = ei{1}.plane{pp}.contexts(condN).rasters.airD.info_metrics.ShannonMI_Zsh;


number_of_columns = size(rasters,2);
column_index = round(timeBefore * effective_sampling_rate);
cis = [];
cis(1,:) = [1 column_index (number_of_columns-column_index+1)];
column_index = cis - 1; column_index(1) = []; column_index(3) = number_of_columns;
cis(2,:) = column_index;
resp = findResponsiveRasters(rasters,cis);
plotRasters_multi({rastersD1,rastersD2,rastersD3,rasters},find(zMIs1'>5),[])
return;
% props = raster_properties(rasters);
% clus = getClusters(props);
n = 0;

n_samples = size(rasters,2);
totalTime = n_samples/effective_sampling_rate;
mdata.xs = linspace(0,totalTime,n_samples);
% % mdata.xs = round(mdata.xs - mdata.xs(column_index),1);
mdata.cis = cis(1,2:3);
plotRasters_simple(rasters,find(resp.ps(:,2)<0.05),mdata)
plotRasters_simple(rastersD,find(resp.ps(:,1)<0.05),[])
plotRasters_multi(rasters,rastersD,find(resp.ps(:,1)<0.05),[])
plotRasters_multi({rastersD1,rastersD2,rastersD3,rasters},find(zMIs1'>5),[])
% plotRasters_simple(rasters,[8 28 19 29 35 47 71 190 197 204 199 186 191 235 282 322 319 344],mdata)

n = 0;