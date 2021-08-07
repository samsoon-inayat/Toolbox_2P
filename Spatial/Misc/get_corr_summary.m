function [adj_PV_corr] = get_corr_summary(all_PV_corr,mask)

N = size(all_PV_corr{1},2);
% mean_corr = cell(N,N);
% mask = ones(size(mean_corr)); 
% mask = triu(mask,order) & ~triu(mask,order+1);
ind = 1;
for rr = 1:N
    for cc = 1:N
%         for_mean_PV_corr = [];
%         for an = 1:length(all_PV_corr)
%             for_mean_PV_corr(:,:,an) = all_PV_corr{an}{rr,cc};
%         end
%         mean_PV_corr{rr,cc} = nanmean(for_mean_PV_corr,3);
        
        % for across adjacent conditions
        if mask(rr,cc)
            for an = 1:length(all_PV_corr)
                this_PV_corr = all_PV_corr{an}{rr,cc};
                PV_corr_diag = diag(this_PV_corr);
                adj_PV_corr{an,ind} = PV_corr_diag;
            end
            ind = ind + 1;
        end
    end
end