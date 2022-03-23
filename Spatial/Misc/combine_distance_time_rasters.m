function out = combine_distance_time_rasters(RsD,RsT,ntrials)

for rr = 1:size(RsD,1)
    for cc = 1:size(RsD,2)
        tRsD = RsD{rr,cc}.sp_rasters1;
        tRsT = RsT{rr,cc}.sp_rasters1;
        szD = size(tRsD,2); szT = size(tRsT,2);
%         tRs = cat(2,tRsD,tRsT);
        mtRsD = squeeze(nanmean(tRsD,1));
        mtRsT = squeeze(nanmean(tRsT,1));
        groups = [ones(1,szD) 2*ones(1,szT)];
        rasters = tRsD;
        p = NaN(size(rasters,3),1);
        resp = logical(zeros(size(rasters,3),1));
        exc = resp; inh = resp;
        parfor ii = 1:size(rasters,3)
            [~,p(ii),~] = ttest2(mtRsD(:,ii),mtRsT(:,ii));
            if p(ii) < 0.05% & hv(ii) == 1
                resp(ii) = 1;
                if nanmean(mtRsD(:,ii)) > nanmean(mtRsT(:,ii))
                    inh(ii) = 1;
                else
                    exc(ii) = 1;
                end
            end
        end
        out.resp{rr,cc} = resp;
        out.exc{rr,cc} = exc;
        out.inh{rr,cc} = inh;
    end
end
n = 0;
propsD = get_props_Rs(RsD,ntrials);
propsT = get_props_Rs(RsT,ntrials);
out.exc_not_good_FR = cell_list_op(out.exc,cell_list_op(propsT.good_FR,[],'not'),'and');
out.inh_not_good_FR = cell_list_op(out.inh,cell_list_op(propsD.good_FR,[],'not'),'and');
out.exc_and_good_FR = cell_list_op(out.exc,propsT.good_FR,'and');
out.inh_and_good_FR = cell_list_op(out.inh,propsD.good_FR,'and');

