function Rs = find_responsive_rasters(Rs,trials)
for rr = 1:size(Rs,1)
    for cc = 1:size(Rs,2)
        R = Rs{rr,cc};
        if strcmp(R.marker_name,'light22T')
            [Rs{rr,cc}.resp.vals,Rs{rr,cc}.resp.cis] = find_resp_time_raster_light(R,trials)
        end
        if strcmp(R.marker_name,'air55T')
            [Rs{rr,cc}.resp.vals,Rs{rr,cc}.resp.cis] = find_resp_time_raster_air(R,trials)
        end
        if strcmp(R.marker_name,'airD') || strcmp(R.marker_name,'beltD')
            [rs,coeffs] = getMRFS_vals(R.gauss_fit_on_mean);
            zMIs = R.info_metrics.ShannonMI_Zsh;
            Rs{rr,cc}.resp.vals = R.iscell' & zMIs > 1.96 & rs > 0.3;
        end
        Rs{rr,cc}.resp.fraction = sum(Rs{rr,cc}.resp.vals)/length(Rs{rr,cc}.resp.vals);
    end
end

    
function [resp,cis] = find_resp_time_raster_air(R,trials)
% SR = R.thorexp.frameRate;
SR = 1/R.bin_width;
markerType = R.marker_name;
timeBefore = str2num(markerType(end-1));
% rasters = R.fromFrames.sp_rasters;
rasters = R.sp_rasters1;
number_of_columns = size(rasters,2);
column_index = round(timeBefore * SR);
cis = [];
cis(1,:) = [1 column_index (number_of_columns-column_index+1)];
column_index = cis - 1; column_index(1) = []; column_index(3) = number_of_columns;
cis(2,:) = column_index;

group = [];
for ii = 1:size(cis,2)
    group = [group ii*ones(1,length(cis(1,ii):cis(2,ii)))];
end
p = NaN(size(rasters,3),1);
p1 = NaN(size(rasters,3),1);
parfor ii = 1:size(rasters,3)
    thisRaster = rasters(trials,:,ii);
    m_thisRaster = nanmean(thisRaster);
%     [p(ii),atab,stats] = anova1(thisRaster,group,'nodisplay');
    [p(ii),atabk,statsk] = kruskalwallis(thisRaster,group,'nodisplay');
    [p1(ii),~,~] = kruskalwallis(m_thisRaster,group,'nodisplay');
end
resp = p < 0.05 & p1 < 0.05;
% resp = p<0.05 & (R.info_metrics.ShannonMI_Zsh > 1.96)';
% % resp = (R.info_metrics.ShannonMI_Zsh > 1.96)';
% [rs,coeffs] = getMRFS_vals(R.gauss_fit_on_mean);
% zMIs = R.info_metrics.ShannonMI_Zsh;
% resp = R.iscell' & zMIs > 1.96 & rs > 0.3;


function [resp,cis] = find_resp_time_raster_light(R,trials)
% SR = R.thorexp.frameRate;
SR = 1/R.bin_width;
markerType = R.marker_name;
timeBefore = str2num(markerType(end-1));
% rasters = R.fromFrames.sp_rasters;
rasters = R.sp_rasters1;
number_of_columns = size(rasters,2);
column_index = round(timeBefore * SR);
cis = [];
cis(1,:) = [1 column_index (number_of_columns-column_index+1)];
column_index = cis - 1; column_index(1) = []; column_index(3) = number_of_columns;
cis(2,:) = column_index;

group = [];
for ii = 1:size(cis,2)
    group = [group ii*ones(1,length(cis(1,ii):cis(2,ii)))];
end
group(group==3) = 2;
p = NaN(size(rasters,3),1);
p1 = NaN(size(rasters,3),1);
parfor ii = 1:size(rasters,3)
    thisRaster = rasters(trials,:,ii);
    m_thisRaster = nanmean(thisRaster);
%     [p(ii),atab,stats] = anova1(thisRaster,group,'nodisplay');
    [p(ii),atabk,statsk] = kruskalwallis(thisRaster,group,'nodisplay');
    [p1(ii),~,~] = kruskalwallis(m_thisRaster,group,'nodisplay');
end
resp = p < 0.05 & p1 < 0.05;
% resp = p<0.05 & (R.info_metrics.ShannonMI_Zsh > 1.96)';
% % resp = (R.info_metrics.ShannonMI_Zsh > 1.96)';
% [rs,coeffs] = getMRFS_vals(R.gauss_fit_on_mean);
% zMIs = R.info_metrics.ShannonMI_Zsh;
% resp = R.iscell' & zMIs > 1.96 & rs > 0.3;
