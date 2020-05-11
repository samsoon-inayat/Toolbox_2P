function SI = spatial_information(ts,spSignal,spj,lj,tj,method)

if strcmp(method,'Dun')
    fi = nanmean(lj./tj);
    f = mean(fi);
    pxi = nanmean(tj)./sum(nanmean(tj));
    SIs = pxi.*(fi./f).*log2(fi./f);
    SI = nansum(SIs);
    return;
end

if strcmp(method,'DunM')
    fi = lj./tj;
    f = nanmean(fi,2);
    pxi = tj./nansum(tj,2);
    SIs = pxi.*(fi./f).*log2(fi./f);
    aSI = nansum(SIs');
    SI = mean(aSI);
    return;
end

if strcmp(method,'MI_formula')
    PS1xj = (lj.*tj)./nansum(lj.*tj,2);
    Pxj = tj./nansum(tj,2);
    PS1 = nansum(lj.*tj,2)./trapz(ts,spSignal);
    SIs = nansum(PS1xj.*Pxj.*log2(PS1xj./PS1),2);
    SI = mode(SIs);
    return;
end

if strcmp(method,'MI_formula1')
    PS1xj = spj./nansum(spj,2);
    Pxj = tj./nansum(tj,2);
    PS1 = nansum(lj.*tj,2)./trapz(ts,spSignal);
    SIs = nansum(PS1xj.*Pxj.*log2(PS1xj./PS1),2);
    SI = mode(SIs);
    return;
end