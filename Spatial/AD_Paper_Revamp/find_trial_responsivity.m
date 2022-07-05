function [presp,all_presp] = find_trial_responsivity(mR)

n = 0;
all_presp = [];
for rr = 1:size(mR,1)
    for cc = 1:size(mR,2)
        tmR = mR{rr,cc};
        for rrr = 1:size(tmR,1)
            for ccc = 1:size(tmR,2)
                ttmR = tmR{rrr,ccc};
                resp = sum(ttmR,2)>0;
                prsp = 100*sum(resp)/length(resp);
                tpresp(rrr,ccc) = prsp;
            end
        end
        presp{rr,cc} = tpresp;
        all_presp = [all_presp tpresp];
    end
end
