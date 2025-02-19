function o = get_props_Rs(Rs,ntrials)

if ~exist('ntrials','var')
    ntrials = 50;
end

% if ~exist('scale','var')
%     scale = NaN;
% end

for rr = 1:size(Rs,1)
    for cc = 1:size(Rs,2)
        if cc == 9
            n = 0;
        end
        R = Rs{rr,cc};
        if size(R.sp_rasters1,2) == 0
            continue;
        end
        nbins = size(R.sp_rasters,2);
        if 50-size(R.speed,2) < 0
          tempS = nanmean(R.speed);
          o.speed{rr,cc} = tempS(1:50);
        else
          o.speed{rr,cc} = padarray(nanmean(R.speed),[0 50-size(R.speed,2)],'post');
        end
        temp_mean = squeeze(nanmean(R.sp_rasters1,[1]));
        o.mean_FR{rr,cc} = (nanmean(temp_mean,1))';
        o.max_FR{rr,cc} = (max(temp_mean,[],1))';
        o.all{rr,cc} = logical(ones(size(R.info_metrics.ShannonMI_Zsh')));
%         if isempty(strfind(R.marker_name,'motion'))
            o.zMI{rr,cc} = R.info_metrics.ShannonMI_Zsh';
            o.MI{rr,cc} = R.info_metrics.ShannonMI';
            try
                o.zMI_MC{rr,cc} = R.info_metrics_MC.ShannonMI_Zsh';
            catch
                o.zMI_MC{rr,cc} = NaN(size(o.zMI{rr,cc}));
            end
            try
                o.MI_MC{rr,cc} = R.info_metrics_MC.ShannonMI';
            catch
                o.MI_MC{rr,cc} = o.zMI_MC{rr,cc};
            end
%         else
%             o.good_zMI_FR{rr,cc} = R.resp.vals & R.resp.FR_based';
%             continue;
%         end
        [rs,MFR,centers,PWs] = get_gauss_fit_parameters(R.gauss_fit_on_mean,R.bin_width);
        o.rs{rr,cc} = rs'; o.MFR{rr,cc} = MFR';  o.PWs{rr,cc} = PWs';
        o.peak_locations{rr,cc} = R.peak_location';
        o.peak_locations_trials{rr,cc} = R.peak_location_trials';
        if strcmp(R.marker_name,'airD')
          centersD = centers';
          o.centersD{rr,cc} = centersD;
          mTime = R.timeStart+(R.timeEnd-R.timeStart)/2;
          firstCol = repmat(mTime(:,1),1,size(mTime,2));
          mTime = mTime - firstCol;
          mTime = mean(mTime);
          for cccc = 1:length(centers)
             ind = find(R.xs - centersD(cccc)>0,1,'first');
             if isempty(ind)
               centers(cccc) = NaN;
             else
               centers(cccc) = mTime(ind);
             end
          end
        end
        if strcmp(R.marker_name,'airID')
          centersD = centers';
          o.centersD{rr,cc} = centersD;
          mTime = R.timeStart+(R.timeEnd-R.timeStart)/2;
          firstCol = repmat(mTime(:,1),1,size(mTime,2));
          mTime = mTime - firstCol;
          mTime = mean(mTime);
          for cccc = 1:length(centers)
             ind = find(R.xs - centersD(cccc)>0,1,'first');
             if isempty(ind)
               centers(cccc) = NaN;
             else
               centers(cccc) = mTime(ind);
             end
          end
        end
        o.centers{rr,cc} = centers';
        
        if strcmp(R.marker_name,'airD')
            bins = 0:50:150;
            bins = 0:75:150;
%             bins = 0:37.5:150;
        end
        if strcmp(R.marker_name,'airIT')
            bins = 0:5:15;
            bins = 0:7.5:15;
%             bins = 0:3.75:15;
        end
        if strcmp(R.marker_name,'airD') || strcmp(R.marker_name,'airIT')
            tpl = o.peak_locations{rr,cc};
            binN = NaN(size(tpl));
            for bii = 1:(length(bins)-1)
                indspeaks = find(tpl >= bins(bii) & tpl < bins(bii+1));
                binN(indspeaks) = bii;
            end
            o.peak_location_bin{rr,cc} = binN;
        end
        
        o.HaFD{rr,cc} = R.fractal_dim.HaFD';
        o.HiFD{rr,cc} = R.fractal_dim.HiFD';
        o.MI_trials_mean{rr,cc} = nanmean(R.MI_trials);
        o.MI_trials{rr,cc} = R.MI_trials;
        xs = R.xs;
%         if ~isempty(strfind(R.marker_name,'D'))
%             p = rs > 0.25 & PWs > xs(2) & PWs < xs(end) & centers >= xs(1)  & centers <= xs(end) & MFR < 10000;
            p = rs > 0.3 & PWs > 1 & PWs < 150 & centers >= 1  & centers <= 150;
            p = PWs < 15;% & centers >= 1  & centers <= 150;
%         end
%         if ~isempty(strfind(R.marker_name,'T'))
%             p = (ones(size(o.zMI{rr,cc})))';
%         end
        o.good_MFR{rr,cc} = MFR' < 10000;
        o.good_Gauss{rr,cc} = p';
        o.good_Gauss_loose{rr,cc} = rs' > 0.25;
        o.good_zMI{rr,cc} = o.zMI{rr,cc} > 1.65;
        o.good_rs{rr,cc} = o.rs{rr,cc} > 0.3;
        if strcmp(R.marker_name,'airIT')
            o.good_PWs{rr,cc} = PWs > 0 & PWs < 15;
            o.good_centers{rr,cc} = centers >=1 & centers <= 15;
        end
        if strcmp(R.marker_name,'airD')
            o.good_PWs{rr,cc} = PWs > 0 & PWs < 150;
            o.good_centers{rr,cc} = centers >=1 & centers <= 150;
        end
        o.nan_zMI{rr,cc} = isnan(o.zMI{rr,cc});
        o.nan_rs{rr,cc} = isnan(o.rs{rr,cc});
        o.good_zMI_MC{rr,cc} = o.zMI_MC{rr,cc} > 1.65;
        [o.good_FR{rr,cc},o.N_Resp_Trials{rr,cc},o.clus_based{rr,cc},o.oc(rr,cc)] = get_FR_based(R.sp_rasters1,ntrials);
        if strcmp(R.marker_name,'airIT')
            [o.good_FR_IT{rr,cc},o.N_Resp_Trials_IT{rr,cc}] = get_FR_based_IT(R.sp_rasters1,ntrials);
            [o.good_FR_IT1{rr,cc},o.N_Resp_Trials_IT1{rr,cc}] = get_FR_based_IT1(R.sp_rasters1,ntrials);
        end
        if strcmp(R.marker_name,'airD')
            [o.good_FR_T{rr,cc},o.N_Resp_Trials_T{rr,cc}] = get_FR_based_T(R.sp_rasters1,ntrials);
            [o.good_FR_T1{rr,cc},o.N_Resp_Trials_T1{rr,cc}] = get_FR_based_T1(R.sp_rasters1,ntrials);
        end
        o.good_FR_and_zMI{rr,cc} = o.good_zMI{rr,cc} & o.good_FR{rr,cc};
        o.good_FR_and_zMI_MC{rr,cc} = o.good_zMI_MC{rr,cc} & o.good_FR{rr,cc};
        temp_z = cell_list_op(o.good_FR(rr,cc),cell_list_op(o.good_zMI(rr,cc),[],'not'),'and');
        o.good_FR_and_notzMI{rr,cc} = temp_z{1};
        
%         o.centers{rr,cc}(~o.good_Gauss{rr,cc}) = NaN;
%         o.PWs{rr,cc}(~o.good_Gauss{rr,cc}) = NaN;
%         o.MFR{rr,cc}(~o.good_Gauss{rr,cc}) = NaN;
        o.trial_scores{rr,cc} = R.resp.trial_scores';
        if isfield(R.resp,'valsA')
            o.valsA{rr,cc} = R.resp.valsA;
        end
        o.valsAD1{rr,cc} = o.good_zMI{rr,cc};
        o.valsAD2{rr,cc} = o.valsAD1{rr,cc} & o.good_MFR{rr,cc};
        o.valsAD3{rr,cc} = o.valsAD2{rr,cc} & o.good_Gauss{rr,cc};
        if size(R.resp.vals,2) == 1
            o.vals{rr,cc} = R.resp.vals;
            o.valsKW{rr,cc} = R.resp.vals;
        else
%             o.vals{rr,cc} = R.resp.vals(:,1);%sum(R.resp.vals,2)>0;
            o.vals{rr,cc} = sum(R.resp.vals,2)>0;
            o.valsKW{rr,cc} = sum(R.resp.valsKW,2)>0;
%             o.vals{rr,cc} = sum(R.resp.vals,2)>floor(size(R.resp.vals,2)/2);
%             tempRespFacVals = R.resp.vals(:,1) | R.resp.vals(:,2);
%             indfac = floor(length(R.resp.fac)/2);
%             if scale == 1
%                 o.vals{rr,cc} = tempRespFacVals & ~R.resp.vals(:,indfac) & ~R.resp.vals(:,end);
%             end
%             if scale == 2
%                 o.vals{rr,cc} = ~tempRespFacVals & R.resp.vals(:,indfac) & ~R.resp.vals(:,end);
%             end
%             if scale == 3
%                 o.vals{rr,cc} = ~tempRespFacVals & R.resp.vals(:,indfac) & R.resp.vals(:,end);
%             end
%             if scale == 4
%                 o.vals{rr,cc} = ~tempRespFacVals & R.resp.vals(:,indfac) & R.resp.vals(:,end);
%             end
        end
        if isfield(R.resp,'valsT')
            if size(R.resp.valsT,2) == 1
                o.valsT{rr,cc} = R.resp.vals;
                o.valsKWT{rr,cc} = R.resp.vals;
            else
                o.valsT{rr,cc} = sum(R.resp.valsT,2)>0;
                o.valsKWT{rr,cc} = sum(R.resp.valsKWT,2)>0;
            end
        end
        temp_cl = cell_list_op(o.good_FR(rr,cc),o.vals(rr,cc),'and');
        temp_c2 = cell_list_op(o.good_FR(rr,cc),cell_list_op(o.vals(rr,cc),[],'not'),'and');
        o.good_FR_and_tuned{rr,cc} = temp_cl{1};
        o.good_FR_and_untuned{rr,cc} = temp_c2{1};
        o.good_HaFD{rr,cc} = o.HaFD{rr,cc} > 1;
        o.good_HiFD{rr,cc} = o.HiFD{rr,cc} > 1;
        o.bad_FR{rr,cc} = ~o.good_FR{rr,cc};
        o.good_FR_and_Gauss{rr,cc} = o.good_FR{rr,cc} & o.good_Gauss{rr,cc};
        temp1 = cell_list_op(o.good_FR(rr,cc),cell_list_op(o.good_Gauss(rr,cc),[],'not'),'and');
        o.good_FR_and_notGauss{rr,cc} = temp1{1};
        
        o.good_FR_and_Gauss_loose{rr,cc} = o.good_FR{rr,cc} & o.good_Gauss_loose{rr,cc};
        temp1 = cell_list_op(o.good_FR(rr,cc),cell_list_op(o.good_Gauss_loose(rr,cc),[],'not'),'and');
        o.good_FR_and_notGauss_loose{rr,cc} = temp1{1};
        
        o.silent_cells{rr,cc} = (o.N_Resp_Trials{rr,cc} == 0);
        if isfield(R.resp,'excinh')
            o.exc{rr,cc} = R.resp.excinh == 1;
            o.inh{rr,cc} = R.resp.excinh == 0;
            temp_cl = cell_list_op(o.good_FR(rr,cc),o.exc(rr,cc),'and');
            temp_c2 = cell_list_op(o.good_FR(rr,cc),o.inh(rr,cc),'and');
            o.good_FR_and_exc{rr,cc} = temp_cl{1};
            o.good_FR_and_inh{rr,cc} = temp_c2{1};
            if strfind(R.marker_name,'motion')
                 temp_m11 = cell_list_op(o.exc(rr,cc),o.inh(rr,cc),'or');
                 o.good_FR{rr,cc} = temp_m11{1};
            end
        end
    end
end
% o.vals_and_good_FR = cell_list_op(o.vals,o.good_FR,'and');
% o.vals_and_good_zMI = cell_list_op(o.vals,o.good_zMI,'and');
% o.not_vals = cell_list_op(o.vals,[],'not'); o.not_good_zMI = cell_list_op(o.good_zMI,[],'not');
% o.vals_and_not_good_zMI = cell_list_op(o.vals,o.not_good_zMI,'and');
% o.not_vals_and_good_zMI = cell_list_op(o.not_vals,o.good_zMI,'and');
% o.not_vals_and_not_good_zMI = cell_list_op(o.not_vals,o.not_good_zMI,'and');
% o.vals_and_valsT = cell_list_op(o.vals,o.valsT,'and');
% o.vals_or_valsT = cell_list_op(o.vals,o.valsT,'or');
% 
% o.valsKW_and_good_FR = cell_list_op(o.valsKW,o.good_FR,'and');
% o.valsKW_and_good_zMI = cell_list_op(o.valsKW,o.good_zMI,'and');
% o.not_valsKW = cell_list_op(o.valsKW,[],'not'); o.not_good_zMI = cell_list_op(o.good_zMI,[],'not');
% o.valsKW_and_not_good_zMI = cell_list_op(o.valsKW,o.not_good_zMI,'and');
% o.not_valsKW_and_good_zMI = cell_list_op(o.not_valsKW,o.good_zMI,'and');
% o.not_valsKW_and_not_good_zMI = cell_list_op(o.not_valsKW,o.not_good_zMI,'and');
% o.valsKW_and_valsKWT = cell_list_op(o.valsKW,o.valsKWT,'and');
% o.valsKW_or_valsKWT = cell_list_op(o.valsKW,o.valsKWT,'or');




function [FR_based,sR,sRp1,oc] = get_FR_based(rasters,ntrials)
rasters = permute(rasters,[2 1 3]);
% sR = sum(squeeze(nanmean(rasters,1))>0.01); % sum over bins and then see in how many trials the sum is greater than 0 to see if the cell responded in multiple trials
sR = sum(squeeze(nansum(rasters,1))>0); % sum over bins and then see in how many trials the sum is greater than 0 to see if the cell responded in multiple trials
sR = sR';
sR = 100*sR./size(rasters,2);
if length(ntrials) == 1
    FR_based = (sR)>=ntrials; % see if cell responded in at least ntrials.
else
    FR_based = (sR)>=ntrials(1) & (sR)<=ntrials(2); % see if cell responded in at least ntrials.
end
% FR_based = FR_based';
trialR = (squeeze(nansum(rasters,1))')>0;
% oc = find_cells_based_on_cluster(trialR);
% sRp = kmeans(trialR,oc);
oc = 2;
rng(1);
sRp = kmeans(trialR,2);
s1 = median(sum(trialR(sRp==1,:),2));
s2 = median(sum(trialR(sRp==2,:),2));
if s1>s2
    sRp1 = sRp == 1;
else
    sRp1 = sRp == 2;
end
% FR_based = sRp1;
n = 0;

function [FR_based,sR,sRp] = get_FR_based_IT(rasters,ntrials)
rasters = rasters(:,1:45,:);
rasters = permute(rasters,[2 1 3]);
% sR = sum(squeeze(nanmean(rasters,1))>0.01); % sum over bins and then see in how many trials the sum is greater than 0 to see if the cell responded in multiple trials
sR = sum(squeeze(nansum(rasters,1))>0); % sum over bins and then see in how many trials the sum is greater than 0 to see if the cell responded in multiple trials
sR = sR';
sR = 100*sR./size(rasters,2);
if length(ntrials) == 1
    FR_based = (sR)>=ntrials; % see if cell responded in at least ntrials.
else
    FR_based = (sR)>=ntrials(1) & (sR)<=ntrials(2); % see if cell responded in at least ntrials.
end
% FR_based = FR_based';

function [FR_based,sR,sRp] = get_FR_based_IT1(rasters,ntrials)
rasters = rasters(:,46:end,:);
rasters = permute(rasters,[2 1 3]);
% sR = sum(squeeze(nanmean(rasters,1))>0.01); % sum over bins and then see in how many trials the sum is greater than 0 to see if the cell responded in multiple trials
sR = sum(squeeze(nansum(rasters,1))>0); % sum over bins and then see in how many trials the sum is greater than 0 to see if the cell responded in multiple trials
sR = sR';
sR = 100*sR./size(rasters,2);
if length(ntrials) == 1
    FR_based = (sR)>=ntrials; % see if cell responded in at least ntrials.
else
    FR_based = (sR)>=ntrials(1) & (sR)<=ntrials(2); % see if cell responded in at least ntrials.
end
% FR_based = FR_based';

function [FR_based,sR,sRp] = get_FR_based_T(rasters,ntrials)
rasters = rasters(:,1:15,:);
rasters = permute(rasters,[2 1 3]);
% sR = sum(squeeze(nanmean(rasters,1))>0.01); % sum over bins and then see in how many trials the sum is greater than 0 to see if the cell responded in multiple trials
sR = sum(squeeze(nansum(rasters,1))>0); % sum over bins and then see in how many trials the sum is greater than 0 to see if the cell responded in multiple trials
sR = sR';
sR = 100*sR./size(rasters,2);
if length(ntrials) == 1
    FR_based = (sR)>=ntrials; % see if cell responded in at least ntrials.
else
    FR_based = (sR)>=ntrials(1) & (sR)<=ntrials(2); % see if cell responded in at least ntrials.
end
% FR_based = FR_based';

function [FR_based,sR,sRp] = get_FR_based_T1(rasters,ntrials)
rasters = rasters(:,16:end,:);
rasters = permute(rasters,[2 1 3]);
% sR = sum(squeeze(nanmean(rasters,1))>0.01); % sum over bins and then see in how many trials the sum is greater than 0 to see if the cell responded in multiple trials
sR = sum(squeeze(nansum(rasters,1))>0); % sum over bins and then see in how many trials the sum is greater than 0 to see if the cell responded in multiple trials
sR = sR';
sR = 100*sR./size(rasters,2);
if length(ntrials) == 1
    FR_based = (sR)>=ntrials; % see if cell responded in at least ntrials.
else
    FR_based = (sR)>=ntrials(1) & (sR)<=ntrials(2); % see if cell responded in at least ntrials.
end
% FR_based = FR_based';

