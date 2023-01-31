function out = find_shifts(allresp_trials,allpeakL_trials,respDT)

for an = 1:5
    for cn = 1:size(allresp_trials,2)
        pLs = allpeakL_trials{an,cn};
        resp = allresp_trials{an,cn};
        if size(respDT,2) == 1
            resp = sum(resp,2)>2 & respDT{an};
        else
            resp = sum(resp,2)>2 & respDT{an,cn};
        end
        pLs = pLs(resp,:);
        [~,inds] = sort(pLs(:,1));
    %     figure(1000);clf;imagesc(pLs(inds,:));colorbar;
%             figure(1000);clf;imagesc(pLs);colorbar;
        100*sum(~isnan(pLs))/size(pLs,1);
        prob_participation = [];
        mshift = []; 
        for ii = 1:size(pLs,1)
            tcell = pLs(ii,:);
            prob_participation(ii,1) = sum(~isnan(tcell))/length(tcell);
            ntcell = tcell(~isnan(tcell));
            mshift(ii,1) = mean(diff(ntcell));
        end
        mshifts(an,cn) = mean(mshift);
        mshiftsB(an,cn) = mean(mshift(mshift<0));
        mshiftsF(an,cn) = mean(mshift(mshift>0));
        fshiftB(an,cn) = sum(mshift<0)/length(mshift);
        fshiftF(an,cn) = sum(mshift>0)/length(mshift);
        fshiftZ(an,cn) = sum(mshift==0)/length(mshift);
    end
end
    
out.mshifts = mshifts;
out.mshiftsB = mshiftsB;
out.mshiftsF = mshiftsF;
out.fshiftsB = fshiftB;
out.fshiftsF = fshiftF;
out.fshiftsZ = fshiftZ;