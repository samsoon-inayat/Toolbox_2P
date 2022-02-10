function Rs = find_responsive_rasters_1(Rs,trialsi)

for rr = 1:size(Rs,1)
    for cc = 1:size(Rs,2)
        R = Rs{rr,cc};
        if isempty(trialsi)
            trials = 1:size(R.sp_rasters,1);
        else
            trials = trialsi;
        end
        if strcmp(R.marker_name,'motionOnsets') || strcmp(R.marker_name,'motionOffsets')
            [Rs{rr,cc}.resp.vals,Rs{rr,cc}.resp.cis,Rs{rr,cc}.resp.excinh] = find_resp_around_event(R,trials,1.5);

        end
        if strcmp(R.marker_name,'light22T') || strcmp(R.marker_name,'tone22T')|| strcmp(R.marker_name,'airOnsets22T') || strcmp(R.marker_name,'airOffsets22T')
            [Rs{rr,cc}.resp.vals,Rs{rr,cc}.resp.cis,Rs{rr,cc}.resp.excinh] = find_resp_around_event(R,trials,2);
        end
        if strcmp(R.marker_name,'airT') || strcmp(R.marker_name,'beltT') || strcmp(R.marker_name,'airIRT') || strcmp(R.marker_name,'airRT')
            zMIs = R.info_metrics.ShannonMI_Zsh';
            if size(R.iscell,2) == 1
                Rs{rr,cc}.resp.vals = R.iscell & zMIs > 1.65;% & R.resp.FR_based;
            else
                Rs{rr,cc}.resp.vals = R.iscell' & zMIs > 1.65;% & R.resp.FR_based;
            end
        end
        if strcmp(R.marker_name,'airIT')
            [Rs{rr,cc}.resp.vals,Rs{rr,cc}.resp.cis] = find_resp_time_raster_intertrial(R,trials);
        end
        if strcmp(R.marker_name,'airD') || strcmp(R.marker_name,'beltD') || strcmp(R.marker_name,'airID') || strcmp(R.marker_name,'airT') || strcmp(R.marker_name,'beltT') || strcmp(R.marker_name,'airIRT') || strcmp(R.marker_name,'airRT') || strcmp(R.marker_name,'airIT')
            [Rs{rr,cc}.resp.vals,Rs{rr,cc}.resp.cis,Rs{rr,cc}.resp.excinh] = find_resp(R,trials);
        end
        Rs{rr,cc}.resp.fraction = sum(Rs{rr,cc}.resp.vals)/length(Rs{rr,cc}.resp.vals);
        Rs{rr,cc}.resp.vals = Rs{rr,cc}.resp.vals;% & Rs{rr,cc}.resp.FR_based;
        n = 0;
    end
end

function [resp,cis,excinh] = find_resp_around_event(R,trials,timeBefore)
% SR = R.thorexp.frameRate;
SR = 1/R.bin_width;
markerType = R.marker_name;
if strfind(markerType,'air')
    n = 0;
end
% timeBefore = str2num(markerType(end-1));
% rasters = R.fromFrames.sp_rasters;
rasters = R.sp_rasters1;
number_of_columns = size(rasters,2);
column_index = round(timeBefore * SR);
cis = [];
cis(1,:) = [1 column_index (number_of_columns-column_index+1)];
column_index = cis - 1; column_index(1) = []; column_index(3) = number_of_columns;
cis(2,:) = column_index;

% resp = (R.info_metrics.ShannonMI_Zsh > 0)';
% return;
group = [];
for ii = 1:size(cis,2)
    group = [group ii*ones(1,length(cis(1,ii):cis(2,ii)))];
end
group(group==3) = 2;
p = NaN(size(rasters,3),1);
p1 = NaN(size(rasters,3),1);
hv = NaN(size(rasters,3),1);
excinh = p;
resp = logical(zeros(size(p)));
parfor ii = 1:size(rasters,3)
    thisRaster = rasters(trials,:,ii);
    m_thisRaster = nanmean(thisRaster);
%     [p(ii),atab,stats] = anova1(thisRaster,group,'nodisplay');
%     [p(ii),~,~] = kruskalwallis(thisRaster,group,'nodisplay');
%     [p(ii),~,~] = kruskalwallis(m_thisRaster,group,'nodisplay');
%     [p(ii),~,~] = kruskalwallis(m_thisRaster,group,'nodisplay');
%     [p(ii),~] = ranksum(m_thisRaster(find(group==1)),m_thisRaster(find(group==2)));
    [~,p(ii),~] = ttest2(m_thisRaster(find(group==1)),m_thisRaster(find(group==2)));
    vert = nansum(thisRaster,2);
    hv(ii) = sum(vert>0) > 5;
%     [~,CRR] = findPopulationVectorPlot(thisRaster,1:10);
%     hv(ii) = findHaFD(CRR,1:size(CRR,1));
    if p(ii) < 0.05% & hv(ii) == 1
        resp(ii,1) = 1;
        if mean(m_thisRaster(group==1)) > mean(m_thisRaster(group==2))
            excinh(ii,1) = 0;
        else
            excinh(ii,1) = 1;
        end
    end
end
% resp = p < 0.05;% & hv;

function [resp,cis,excinh] = find_resp(R,trials)
% SR = R.thorexp.frameRate;
SR = 1/R.bin_width;
markerType = R.marker_name;
if strfind(markerType,'air')
    n = 0;
end
% timeBefore = str2num(markerType(end-1));
% rasters = R.fromFrames.sp_rasters;
rasters = R.sp_rasters1;
number_of_columns = size(rasters,2);
if mod(number_of_columns,2)
    
end
num_groups = 10;
col_per_group = floor(number_of_columns/num_groups);
num_groups = ceil(number_of_columns/col_per_group);
group = [];
for ii = 1:num_groups
    group = [group ii*ones(1,col_per_group)];
end
group = group(1:number_of_columns);
p = NaN(size(rasters,3),1);
p1 = NaN(size(rasters,3),1);
hv = NaN(size(rasters,3),1);
excinh = p;
resp = logical(zeros(size(p)));
for ii = 1:size(rasters,3)
    thisRaster = rasters(trials,:,ii);
    [within,dvn,xlabels] = make_within_table({'Col'},[number_of_columns]);
    dataT = make_between_table({thisRaster},dvn);
    ra = RMA(dataT,within);
    try
        p(ii) = ra.ranova.pValue_sel(2);
    catch
        p(ii) = 1;
    end
    

%     m_thisRaster = nanmean(thisRaster);
%     [p(ii),atab,stats] = anova1(thisRaster,group,'nodisplay');
%     [c,m,h,gnames] = multcompare(stats);
%     [p(ii),atab,stats] = anova1(thisRaster,1:number_of_columns,'nodisplay');
%     [c,m,h,gnames] = multcompare(stats);
% %     [p(ii),~,~] = kruskalwallis(thisRaster,group,'nodisplay');
% %     [p(ii),~,~] = kruskalwallis(m_thisRaster,group,'nodisplay');
% %     [p(ii),~,~] = kruskalwallis(m_thisRaster,group,'nodisplay');
% %     [p(ii),~] = ranksum(m_thisRaster(find(group==1)),m_thisRaster(find(group==2)));
%     [~,p(ii),~] = ttest2(m_thisRaster(find(group==1)),m_thisRaster(find(group==2)));
%     vert = nansum(thisRaster,2);
%     hv(ii) = sum(vert>0) > 5;
%     [~,CRR] = findPopulationVectorPlot(thisRaster,1:10);
%     hv(ii) = findHaFD(CRR,1:size(CRR,1));
    if p(ii) < 0.05% & hv(ii) == 1
        resp(ii,1) = 1;
%         if mean(m_thisRaster(group==1)) > mean(m_thisRaster(group==2))
%             excinh(ii,1) = 0;
%         else
%             excinh(ii,1) = 1;
%         end
    end
end
% resp = p < 0.05;% & hv;

