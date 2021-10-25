function o = get_props_Rs(Rs,ntrials)

if ~exist('ntrials','var')
    ntrials = 5;
end

for rr = 1:size(Rs,1)
    for cc = 1:size(Rs,2)
        R = Rs{rr,cc};
        o.all{rr,cc} = logical(ones(size(R.info_metrics.ShannonMI_Zsh')));
%         if isempty(strfind(R.marker_name,'motion'))
            o.zMI{rr,cc} = R.info_metrics.ShannonMI_Zsh';
            o.MI{rr,cc} = R.info_metrics.ShannonMI';
%         else
%             o.good_zMI_FR{rr,cc} = R.resp.vals & R.resp.FR_based';
%             continue;
%         end
        [rs,MFR,centers,PWs] = get_gauss_fit_parameters(R.gauss_fit_on_mean,R.bin_width);
        o.rs{rr,cc} = rs'; o.MFR{rr,cc} = MFR'; o.centers{rr,cc} = centers'; o.PWs{rr,cc} = PWs';
        o.peak_locations{rr,cc} = R.peak_location';
        o.peak_locations_trials{rr,cc} = R.peak_location_trials';
        
        o.HaFD{rr,cc} = R.fractal_dim.HaFD';
        o.HiFD{rr,cc} = R.fractal_dim.HiFD';
        
        xs = R.xs;
%         if ~isempty(strfind(R.marker_name,'D'))
            p = rs > 0.25 & PWs > xs(2) & PWs < xs(end) & centers > xs(1)  & centers < xs(end) & MFR < 10000;
%         end
%         if ~isempty(strfind(R.marker_name,'T'))
%             p = (ones(size(o.zMI{rr,cc})))';
%         end
        o.good_Gauss{rr,cc} = p';
        o.good_zMI{rr,cc} = o.zMI{rr,cc} > 1.65;
        [o.good_FR{rr,cc},o.N_Resp_Trials{rr,cc}] = get_FR_based(R.sp_rasters1,ntrials);
        o.good_FR_and_zMI{rr,cc} = o.good_zMI{rr,cc} & o.good_FR{rr,cc};
        o.centers{rr,cc}(~o.good_Gauss{rr,cc}) = NaN;
        o.PWs{rr,cc}(~o.good_Gauss{rr,cc}) = NaN;
        o.MFR{rr,cc}(~o.good_Gauss{rr,cc}) = NaN;
        o.trial_scores{rr,cc} = R.resp.trial_scores';
        o.vals{rr,cc} = R.resp.vals;
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
        o.silent_cells{rr,cc} = (o.N_Resp_Trials{rr,cc} == 0);
        if isfield(R.resp,'excinh')
            o.exc{rr,cc} = R.resp.excinh == 1;
            o.inh{rr,cc} = R.resp.excinh == 0;
            temp_cl = cell_list_op(o.good_FR(rr,cc),o.exc(rr,cc),'and');
            temp_c2 = cell_list_op(o.good_FR(rr,cc),o.inh(rr,cc),'and');
            o.good_FR_and_exc{rr,cc} = temp_cl{1};
            o.good_FR_and_inh{rr,cc} = temp_c2{1};
        end
    end
end

function [FR_based,sR,sRp] = get_FR_based(rasters,ntrials)
rasters = permute(rasters,[2 1 3]);
% sR = sum(squeeze(nanmean(rasters,1))>0.01); % sum over bins and then see in how many trials the sum is greater than 0 to see if the cell responded in multiple trials
sR = sum(squeeze(nansum(rasters,1))>0); % sum over bins and then see in how many trials the sum is greater than 0 to see if the cell responded in multiple trials
sR = sR';
sR = 100*sR./size(rasters,2);
FR_based = (sR)>=ntrials; % see if cell responded in at least ntrials.
% FR_based = FR_based';

