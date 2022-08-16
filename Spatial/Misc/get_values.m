function [all_vals,all_vals_N] = get_values(Rs,number_of_bins,var,sel_pop)

all_vals = [];
all_vals_N = [];

for cc = 1:size(Rs,2)
    these_vals = NaN(size(Rs,1),number_of_bins);
    these_vals_N = NaN(size(Rs,1),number_of_bins);
    for rr = 1:size(Rs,1)
        if cc == 4
            n = 0;
        end
        R = Rs{rr,cc};
        mbl = mean(R.beltLength);
        binSize = mbl/number_of_bins;
        resp = sel_pop{rr,cc};
        [rs,as,bs,cs] = get_gauss_fit_parameters(R.gauss_fit_on_mean,R.bin_width);
        if strcmp(var,'all_zMIs')
            vals = R.info_metrics.ShannonMI_Zsh';
        end
        if strcmp(var,'all_fFR')
            vals = as';
%             resp(vals>5000) = 0;
        end
        if strcmp(var,'all_fwidths')
            vals = cs';
        end
        if strcmp(var,'all_frs')
            vals = rs';
        end
        vals = vals(resp);
        if number_of_bins > 1
        bs = bs(resp);
        bins = 0:binSize:mbl;
        bins = [0 50 100 150];
        [N,E,Bi] = histcounts(bs,bins);
        mean_sv_Vals = [];
        for bb = 1:length(N)
            mean_sv_Vals(bb) = nanmean(vals(Bi == bb));
        end
        these_vals(rr,:) = mean_sv_Vals;
        these_vals_N(rr,:) = 100*(N/sum(N));
        these_vals_N(rr,:) = 100*(N/length(resp));
        else
            these_vals(rr,:) = nanmean(vals);
            these_vals_N(rr,:) = 100*length(vals)/length(resp);
        end
    end
    all_vals = [all_vals these_vals];
    all_vals_N = [all_vals_N these_vals_N];
end
