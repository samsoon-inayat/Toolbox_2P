function [all_OI,mOI,semOI,all_OI_mat] = get_overlap_index(resp_valsi)


% resp_valsCi = resp_valsC;
for rr = 1:size(resp_valsi,1)
    ccs = [];
    for cc = 1:size(resp_valsi,2)
        ccs(:,cc) = resp_valsi{rr,cc};
    end
    resp_vals{rr} = ccs;
end


for ii = 1:length(resp_vals)
    ccs = resp_vals{ii};
    OI = NaN(size(ccs,2),size(ccs,2));
    mask = triu(ones(size(ccs,2)),1);
    for rr = 1:size(ccs,2)
        ccs1 = ccs(:,rr);
        nan1 = isnan(ccs1);
        for cc = 1:size(ccs,2)
%             if rr == cc
%                 continue;
%             end
%             if mask(rr,cc)
                ccs2 = ccs(:,cc);
                nan2 = isnan(ccs2);
                inds_nans = nan1 | nan2;
                [ia1,ib1] = sort(ccs1(~inds_nans));
                [ia2,ib2] = sort(ccs2(~inds_nans));
                sh = [];
                for dd = 1:length(ib1)
                    ind2 = find(ib2 == ib1(dd));
                    sh(dd) = 100*(ind2 - dd)/size(ccs1,1);
                end
                OI(rr,cc) = (mean(abs(sh)));
%             end
        end
    end
    all_OI{ii} = OI;
    all_OI_mat(:,:,ii) = OI;
end

mOI = NaN(size(all_OI{1},1),size(all_OI{1},2),length(all_OI));

for ii = 1:length(all_OI)
    mOI(:,:,ii) = all_OI{ii};
end
semOI = std(mOI,[],3)/sqrt(length(all_OI));
mOI = mean(mOI,3);
