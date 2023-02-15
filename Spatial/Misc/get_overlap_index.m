function [all_OI,mOI,semOI,all_OI_mat,p_vals,h,all_CI,mCI,semCI,all_CI_mat,uni] = get_overlap_index(resp_valsi,thr,alpha,inv)


if ~exist('inv','var')
    inv = 0;
end

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
        for cc = 1:size(ccs,2)
%             if rr == cc
%                 continue;
%             end
%             if mask(rr,cc)
                ccs2 = ccs(:,cc);
                if inv == 1
                    ccs2 = ~ccs2;
                end
                shared = ccs1 & ccs2;
                OI(rr,cc) = sum(shared)/(sum(ccs1)+sum(ccs2)-sum(shared));
                if OI(rr,cc) < 0
                    n = 0;
                end
                CI(rr,cc) = 100*sum(shared)/length(ccs1);
                n = 0;
                tuni(rr,cc) = 100*sum(ccs1 & (~ccs2))/length(ccs1);
%             end
        end
    end
    all_OI{ii} = 1-OI;
    all_OI_mat(:,:,ii) = OI;
    all_CI{ii} = CI;
    all_CI_mat(:,:,ii) = CI;
    uni(:,:,ii) = tuni;
end

mOI = NaN(size(all_OI{1},1),size(all_OI{1},2),length(all_OI));

for ii = 1:length(all_OI)
    mOI(:,:,ii) = all_OI{ii};
end
semOI = std(mOI,[],3)/sqrt(length(all_OI));
mOI = mean(mOI,3);

mCI = NaN(size(mOI));

for ii = 1:length(all_CI)
    mCI(:,:,ii) = all_CI{ii};
end
semCI = std(mCI,[],3)/sqrt(length(all_CI));
mCI = mean(mCI,3);


h = NaN(size(all_OI_mat(:,:,1)));

for rr = 1:size(all_OI_mat,1)
    for cc = 1:size(all_OI_mat,2)
        if isnan(all_OI_mat(rr,cc,1))
            continue;
        end
        [h(rr,cc),p_vals(rr,cc)] = ttest(squeeze(all_OI_mat(rr,cc,:)),thr,'Alpha',alpha);
    end
end
