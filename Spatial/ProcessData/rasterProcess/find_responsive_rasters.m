function Rs = find_responsive_rasters(Rs,trialsi)

for rr = 1:size(Rs,1)
    for cc = 1:size(Rs,2)
        R = Rs{rr,cc};
        rasters = permute(R.sp_rasters1,[2 1 3]);
        sR = sum(squeeze(nanmean(rasters,1))>0); % sum over bins and then see in how many trials the sum is greater than 0 to see if the cell responded in multiple trials
        R.resp.FR_based36 = (sR)>=3 & (sR)<=6; % see if cell responded in at least half trials.
        R.resp.FR_based710 = (sR)>=7 & (sR)<=10; % see if cell responded in at least half trials.
        total_trials = size(rasters,2);
        R.resp.trial_scores = sR/total_trials;
        R.resp.trial_scores_percent = 100*sR/total_trials;
        Rs{rr,cc} = R;
    end
end

for rr = 1:size(Rs,1)
    for cc = 1:size(Rs,2)
        R = Rs{rr,cc};
        if isempty(trialsi)
            trials = 1:size(R.sp_rasters,1);
        else
            trials = trialsi;
        end
        if strcmp(R.marker_name,'motionOnsets') || strcmp(R.marker_name,'motionOffsets')
            [Rs{rr,cc}.resp.vals,Rs{rr,cc}.resp.cis,Rs{rr,cc}.resp.excinh] = find_resp_motionOnset_raster(R,trials);
            if 0
                Rs{rr,cc}.sp_rasters1 = Rs{rr,cc}.fromFrames.sp_rasters;
                Rs{rr,cc}.dist = Rs{rr,cc}.fromFrames.dist;
                Rs{rr,cc}.speed = Rs{rr,cc}.fromFrames.speed;
                Rs{rr,cc}.duration = Rs{rr,cc}.fromFrames.duration;
                Rs{rr,cc}.xs = Rs{rr,cc}.fromFrames.duration(1,:);
                Rs{rr,cc}.bin_width = 1/Rs{rr,cc}.thorexp.frameRate;
            end
            n = 0;
        end
        if strcmp(R.marker_name,'light22T') || strcmp(R.marker_name,'tone22T')|| strcmp(R.marker_name,'airOnsets22T') || strcmp(R.marker_name,'airOffsets22T')
            [Rs{rr,cc}.resp.vals,Rs{rr,cc}.resp.cis,Rs{rr,cc}.resp.excinh] = find_resp_time_raster_light(R,trials);
%             [Rs{rr,cc}.resp.vals,Rs{rr,cc}.resp.cis] = find_resp_time_raster_light_fractal(R,trials);
        end
        if strcmp(R.marker_name,'air55T') ||  strcmp(R.marker_name,'air77T') || strcmp(R.marker_name,'air44T')
            [Rs{rr,cc}.resp.vals,Rs{rr,cc}.resp.cis,Rs{rr,cc}.resp.excinh] = find_resp_time_raster_air(R,trials);
        end
        if strcmp(R.marker_name,'air33T') 
            [Rs{rr,cc}.resp.vals,Rs{rr,cc}.resp.cis] = find_resp_time_raster_air(R,trials);
        end
        if strcmp(R.marker_name,'airT') || strcmp(R.marker_name,'beltT') || strcmp(R.marker_name,'airIRT') || strcmp(R.marker_name,'airRT')
            zMIs = R.info_metrics.ShannonMI_Zsh';
            if size(R.iscell,2) == 1
                Rs{rr,cc}.resp.vals = R.iscell & zMIs > 1.65;% & R.resp.FR_based;
            else
                Rs{rr,cc}.resp.vals = R.iscell' & zMIs > 1.65;% & R.resp.FR_based;
            end
            [Rs{rr,cc}.resp.vals,Rs{rr,cc}.resp.fac] = find_resp_no_before_after_anova(R,1:10);
            Rs{rr,cc}.resp.valsC = sum(Rs{rr,cc}.resp.vals,2)>0;
        end
        if strcmp(R.marker_name,'airIT') || strcmp(R.marker_name,'airOnsets55T') || strcmp(R.marker_name,'airOffsets55T')
%             [Rs{rr,cc}.resp.vals,Rs{rr,cc}.resp.cis] = find_resp_time_raster_intertrial(R,trials);
            [rs,MFR,centers,PWs] = get_gauss_fit_parameters(R.gauss_fit_on_mean,R.bin_width);
            zMIs = R.info_metrics.ShannonMI_Zsh;
%             Rs{rr,cc}.resp.vals = R.iscell' & zMIs > 1.96 & rs > 0.25 & PWs < 150 & centers > 0 & centers < 150;
%             Rs{rr,cc}.resp.vals = R.iscell' & zMIs > 1.65 & rs > 0.25 & PWs > 1 & PWs < 15 & centers > 0 & centers < 15 & MFR < 10000;
            Rs{rr,cc}.resp.vals = R.iscell' & zMIs > 1.65 & rs > 0.25;
            Rs{rr,cc}.resp.vals = Rs{rr,cc}.resp.vals';
            [Rs{rr,cc}.resp.vals,Rs{rr,cc}.resp.fac] = find_resp_no_before_after_anova(R,1:10);
            Rs{rr,cc}.resp.valsC = sum(Rs{rr,cc}.resp.vals,2)>0;
        end
        if strcmp(R.marker_name,'airD') || strcmp(R.marker_name,'beltD') || strcmp(R.marker_name,'airID')
%             [rs1,coeffs] = getMRFS_vals(R.gauss_fit_on_mean);
            [rs,MFR,centers,PWs] = get_gauss_fit_parameters(R.gauss_fit_on_mean,R.bin_width);
            zMIs = R.info_metrics.ShannonMI_Zsh;
%             Rs{rr,cc}.resp.vals = R.iscell' & zMIs > 1.96 & rs > 0.25 & PWs < 150 & centers > 0 & centers < 150;
            Rs{rr,cc}.resp.vals = R.iscell' & zMIs > 1.65 & rs > 0.25 & PWs > 1 & PWs < 150 & centers > 0 & centers < 150 & MFR < 10000;
            Rs{rr,cc}.resp.vals = Rs{rr,cc}.resp.vals';
            [Rs{rr,cc}.resp.vals,Rs{rr,cc}.resp.fac] = find_resp_no_before_after_anova(R,1:10);
            Rs{rr,cc}.resp.valsC = sum(Rs{rr,cc}.resp.vals,2)>0;
        end
        Rs{rr,cc}.resp.fraction = sum(Rs{rr,cc}.resp.vals)/length(Rs{rr,cc}.resp.vals);
%         Rs{rr,cc}.resp.vals = Rs{rr,cc}.resp.vals;% & Rs{rr,cc}.resp.FR_based;
        n = 0;
    end
end

function [resp,cis,excinh] = find_resp_motionOnset_raster(R,trials)
% SR = R.thorexp.frameRate;
% rasters = R.fromFrames.sp_rasters;
SR = 1/R.bin_width;
rasters = R.sp_rasters1;
markerType = R.marker_name;
timeBefore = 1.5;

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
resp = logical(zeros(size(p)));
excinh = p;
for ii = 1:size(rasters,3)
    thisRaster = rasters(trials,:,ii);
    m_thisRaster = nanmean(thisRaster);
%     [p(ii),atab,stats] = anova1(thisRaster,group,'nodisplay');
%     [p(ii),~,~] = kruskalwallis(thisRaster,group,'nodisplay');
%     [p(ii),~,~] = kruskalwallis(m_thisRaster,group,'nodisplay');
%     [p(ii),~] = ranksum(m_thisRaster(find(group==1)),m_thisRaster(find(group==2)));
    [~,p(ii),~] = ttest2(m_thisRaster(find(group==1)),m_thisRaster(find(group==2)));
    vert = nansum(thisRaster,2);
    hv(ii) = sum(vert>0) > (round(length(trials)/2)-1);
%     [~,CRR] = findPopulationVectorPlot(thisRaster,1:10);
%     hv(ii) = findHaFD(CRR,1:size(CRR,1));
    if p(ii) < 0.05% & hv(ii) == 1
        resp(ii) = 1;
        if mean(m_thisRaster(group==1)) > mean(m_thisRaster(group==2))
            excinh(ii) = 0;
        else
            excinh(ii) = 1;
        end
    end
end



function [resp,cis,excinh] = find_resp_time_raster_intertrial(R,trials)
% SR = R.thorexp.frameRate;
SR = 1/R.bin_width;
% markerType = R.marker_name;
timeBefore = 7.5;
% rasters = R.fromFrames.sp_rasters;
rasters = R.sp_rasters1;
number_of_columns = size(rasters,2);
column_index = round(timeBefore * SR);
cis = [];
cis(1,:) = [1 column_index (number_of_columns-column_index+1)];
column_index = cis - 1; column_index(1) = []; column_index(3) = number_of_columns;
cis(2,:) = column_index;

% resp = (R.info_metrics.ShannonMI_Zsh > 1.96)';
% return;

group = [];
for ii = 1:size(cis,2)
    group = [group ii*ones(1,length(cis(1,ii):cis(2,ii)))];
end
group(group == 3) = 2;
p = NaN(size(rasters,3),1);
hv = p;
resp = logical(zeros(size(p)));
excinh = p;
parfor ii = 1:size(rasters,3)
    thisRaster = rasters(trials,:,ii);
    m_thisRaster = nanmean(thisRaster);
    [p(ii),atabk,statsk] = kruskalwallis(m_thisRaster,group,'nodisplay');
    vert = nanmean(thisRaster,2)>0.1;
%     vert = nansum(thisRaster,2);
%     hv(ii) = sum(vert>0) >= 5;
%     hv(ii) = sum(vert) >= 5;
    if p(ii) < 0.05% & hv(ii) == 1
        resp(ii,1) = 1;
        if mean(m_thisRaster(group==1)) > mean(m_thisRaster(group==2)) && mean(m_thisRaster(group==3)) > mean(m_thisRaster(group==2))
            excinh(ii) = 0;
        else
            excinh(ii) = 1;
        end
    end
end



function [resp,cis,excinh] = find_resp_time_raster_air(R,trials)
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

% resp = (R.info_metrics.ShannonMI_Zsh > 1.96)';
% return;

group = [];
for ii = 1:size(cis,2)
    group = [group ii*ones(1,length(cis(1,ii):cis(2,ii)))];
end
p = NaN(size(rasters,3),1);
hv = p;
resp = logical(zeros(size(p)));
excinh = p;
parfor ii = 1:size(rasters,3)
    thisRaster = rasters(trials,:,ii);
    m_thisRaster = nanmean(thisRaster);
%     [p(ii),atabk,statsk] = kruskalwallis(m_thisRaster,group,'nodisplay');
    [p(ii),atabk,statsk] = anova1(m_thisRaster,group,'nodisplay');
    vert = nanmean(thisRaster,2)>0.1;
%     vert = nansum(thisRaster,2);
%     hv(ii) = sum(vert>0) >= 5;
    hv(ii) = sum(vert) >= 5;
    if p(ii) < 0.05% & hv(ii) == 1
        resp(ii,1) = 1;
        if mean(m_thisRaster(group==1)) > mean(m_thisRaster(group==2)) && mean(m_thisRaster(group==3)) > mean(m_thisRaster(group==2))
            excinh(ii) = 0;
        else
            excinh(ii) = 1;
        end
    end
end


% resp = p<0.05 & (R.info_metrics.ShannonMI_Zsh > 1.96)';
% % resp = (R.info_metrics.ShannonMI_Zsh > 1.96)';
% [rs,coeffs] = getMRFS_vals(R.gauss_fit_on_mean);
% zMIs = R.info_metrics.ShannonMI_Zsh;
% resp = zMIs > 1.96 & rs > 0.3;


function [resp,cis,excinh] = find_resp_time_raster_light(R,trials)
% SR = R.thorexp.frameRate;
SR = 1/R.bin_width;
markerType = R.marker_name;
if strfind(markerType,'air')
    n = 0;
end
timeBefore = str2num(markerType(end-1));
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


function [resp,cis] = find_resp_time_raster_light_fractal(R,trials)
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
[rs,as,bs,PWs] = get_gauss_fit_parameters(R.gauss_fit_on_mean,R.bin_width);
RS= rs'; SI = R.info_metrics.ShannonMI_Zsh'; HiFD = R.fractal_dim.HiFD';
HaFD = R.fractal_dim.HaFD';
propMatrix = [SI RS HiFD HaFD];
rng(3,'twister');
num_clus = 50;
[clus_inds,clus_centers] = kmeans(propMatrix,num_clus);
for ii = 1:num_clus
    cluster_median_SI(ii) = median(SI(clus_inds == ii));
    cluster_mean_SI(ii) = mean(SI(clus_inds == ii));
    cluster_min_SI(ii) = min(SI(clus_inds == ii));
    cluster_min_HaFD(ii) = min(HaFD(clus_inds == ii));
    cluster_min_HiFD(ii) = min(HiFD(clus_inds == ii));
    cluster_min_Rs(ii) = min(RS(clus_inds == ii));
end
clus_More = find(cluster_min_Rs > 0.3 & cluster_min_SI > 0);
n_clus_inds = zeros(size(clus_inds));
for ii = 1:length(clus_More)
    n_clus_inds = n_clus_inds | clus_inds == clus_More(ii);
end
resp = n_clus_inds;% & rs' > 0.3;



function [resp,fac] = find_resp_no_before_after_anova(R,trials)
% SR = R.thorexp.frameRate;
file_name = fullfile(R.pd_folder,sprintf('responsive_cells_pyramidal_anova_%s.mat',R.context_info));
if exist(file_name,'file')
    te = load(file_name);
    resp = te.resp;
    resp = te.p<0.05;
    fac = te.fac;
    return;
end
SR = 1/R.bin_width;
markerType = R.marker_name;
timeBefore = str2num(markerType(end-1));
% rasters = R.fromFrames.sp_rasters;
rasters = R.sp_rasters1;
number_of_columns = size(rasters,2);
half_num_col = floor(number_of_columns/2);
fac = 1:1:half_num_col;
for ii = 1:length(fac)
    os = ones(1,fac(ii));
    num_groups = floor(number_of_columns/fac(ii));
    tgroup = [];
    for jj = 1:num_groups
        tgroup = [tgroup jj*os];
    end
    diffL = number_of_columns - length(tgroup);
    if diffL > 0
        tgroup = [tgroup ones(1,diffL)*(tgroup(end)+1)];
    end
    groups{ii} = tgroup;
end
p = NaN(size(rasters,3),length(fac));
resp = logical(zeros(size(p)));
pG = p;
for ni = 1:length(fac)
    group = groups{ni};
    parfor ii = 1:size(rasters,3)
        thisRaster = rasters(trials,:,ii);
%         thisRasterG = make_gauss_fit_raster(R,ii);
        p(ii,ni) = anova1(thisRaster,group,'nodisplay');
%         pG(ii,ni) = anova1(thisRasterG,group,'nodisplay');
    end
end
resp = sum(p<0.05,2)>0;
save(file_name,'resp','p','fac');
% respG = sum(pG<0.05,2)>0;



function [resp,cis,excinh] = find_resp_no_before_after(R,trials)
% SR = R.thorexp.frameRate;
SR = 1/R.bin_width;
markerType = R.marker_name;
timeBefore = str2num(markerType(end-1));
% rasters = R.fromFrames.sp_rasters;
rasters = R.sp_rasters1;
number_of_columns = size(rasters,2);
column_index = round(timeBefore * SR);
cis = [];
num_in_group = ceil(number_of_columns/10);
gs = zeros(1,number_of_columns);
gs(1:num_in_group) = 1;
gs = logical(gs);
for jj = 1:number_of_columns
    group = circshift(gs,jj-1);
    p = NaN(size(rasters,3),1);
    resp = logical(zeros(size(p)));
    excinh = p;
    parfor ii = 1:size(rasters,3)
        thisRaster = rasters(trials,:,ii);
        m_thisRaster = nanmean(thisRaster);
    %     [p(ii),atabk,statsk] = kruskalwallis(m_thisRaster,group,'nodisplay');
    %     [p(ii),atabk,statsk] = anova1(m_thisRaster,group,'nodisplay');
        [~,p(ii),~] = ttest2(m_thisRaster(find(group==1)),m_thisRaster(find(group==0)));
    %     vert = nanmean(thisRaster,2)>0.1;
    %     vert = nansum(thisRaster,2);
    %     hv(ii) = sum(vert>0) >= 5;
    %     hv(ii) = sum(vert) >= 5;
        if p(ii) < 0.05% & hv(ii) == 1
            resp(ii,1) = 1;
            if mean(m_thisRaster(group==1)) > mean(m_thisRaster(group==0))
                excinh(ii) = 0;
            else
                excinh(ii) = 1;
            end
        end
    end
    a_p(:,jj) = p;
    a_excinh(:,jj) = excinh;
    a_resp(:,jj) = resp;
end
resp = sum(a_resp,2)>2;
excinh = sum(a_excinh,2);

