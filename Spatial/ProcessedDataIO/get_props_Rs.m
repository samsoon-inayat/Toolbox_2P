function o = get_props_Rs(Rs,ntrials)

if ~exist('ntrials','var')
    ntrials = 5;
end

for rr = 1:size(Rs,1)
    for cc = 1:size(Rs,2)
        R = Rs{rr,cc};
%         if isempty(strfind(R.marker_name,'motion'))
            o.zMI{rr,cc} = R.info_metrics.ShannonMI_Zsh';
%         else
%             o.good_zMI_FR{rr,cc} = R.resp.vals & R.resp.FR_based';
%             continue;
%         end
        [rs,MFR,centers,PWs] = get_gauss_fit_parameters(R.gauss_fit_on_mean,R.bin_width);
        o.rs{rr,cc} = rs'; o.MFR{rr,cc} = MFR'; o.centers{rr,cc} = centers'; o.PWs{rr,cc} = PWs';
        o.peak_locations{rr,cc} = R.peak_location';
        
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
        o.good_zMI{rr,cc} = o.zMI{rr,cc} > 1.65
        o.good_FR{rr,cc} = get_FR_based(R.sp_rasters1,ntrials);
        o.good_zMI_FR{rr,cc} = o.good_zMI{rr,cc} & o.good_FR{rr,cc};
        o.centers{rr,cc}(~o.good_Gauss{rr,cc}) = NaN;
        o.PWs{rr,cc}(~o.good_Gauss{rr,cc}) = NaN;
        o.MFR{rr,cc}(~o.good_Gauss{rr,cc}) = NaN;
    end
end

function FR_based = get_FR_based(rasters,ntrials)
rasters = permute(rasters,[2 1 3]);
sR = sum(squeeze(nansum(rasters,1))>0); % sum over bins and then see in how many trials the sum is greater than 0 to see if the cell responded in multiple trials
FR_based = (sR)>ntrials; % see if cell responded in at least half trials.
FR_based = FR_based';

