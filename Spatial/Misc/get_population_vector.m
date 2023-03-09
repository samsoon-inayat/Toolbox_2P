function [pv,pvc,pvct,outg] = get_population_vector(mR,resp)


common = 0;
outg = [];
if ~iscell(resp)
    resp = repmat({resp},1,size(mR,2));
    common = 1;
else
    if size(resp,2) < size(mR,2)
        resp = repmat(resp,1,size(mR,2));
        common = 1;
    end
end

rr = 1;
for cc = 1:size(mR,2)
    tmR = mR{rr,cc};
    tpv = tmR(resp{rr,cc},:);
    pv{rr,cc} = tpv;
    pvc{rr,cc} = corr(tpv);
    pvct{rr,cc} = corr(tpv');
end

if common
    gmR = [];
    for cc = 1:size(mR,2)
        tmR = pv{rr,cc};
        gmR = [gmR tmR];
    end
    outg.pv = gmR;
    outg.pvc = corr(gmR);
    outg.pvct = corr(gmR');
end
