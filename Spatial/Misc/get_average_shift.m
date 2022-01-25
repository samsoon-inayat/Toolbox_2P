function [all_OI,mOI,semOI,all_OI_mat,mOIM,semOIM,OIM_mat] = get_overlap_index(resp_valsi)


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
    mmOI(:,:,ii) = all_OI{ii};
end
semOI = std(mOI,[],3)/sqrt(length(all_OI));
mOI = mean(mOI,3);

sinds = 1:10:size(mmOI,1);
einds = 10:10:size(mmOI,1);

for an = 1:size(mmOI,3)
    tmmOI = mmOI(:,:,an);
    mmmOI = [];
    for ii = 1:length(sinds)
        for jj = 1:length(sinds)
            indsr = sinds(ii):einds(ii);
            indsc = sinds(jj):einds(jj);
            mmmOI(ii,jj) = mean(tmmOI(indsr,indsc),'A');
        end
    end
    mmmOIan(:,:,an) = mmmOI;
end

semOIM = std(mmmOIan,[],3)/sqrt(length(all_OI));
mOIM = mean(mmmOIan,3);
OIM_mat = mmmOIan;

n = 0;



