function [rasterF,mdl,rs] = do_gauss_fit(xs,raster,statsetfitnlm,orderdc)
order = orderdc(1);
dc = orderdc(2);
rasterF = NaN(size(raster));
rs = NaN(size(raster,1),1+order*3);
for ii = 1:size(raster,1)
    mSig = (raster(ii,:));
    [max_mSig,loc_max_mSig] = max(mSig);
    [pks,locs,w,p] = findpeaks(mSig);
    inds = p > nanmean(p);
    if sum(inds) > 1
        pks = pks(inds);locs = locs(inds);w = w(inds);
    end
    Rs = NaN(1,length(pks));
    for jj = 1:length(pks)
        thisC = [pks(jj),locs(jj),w(jj)];
        [gmdl,x0] = buildGaussianModel(order,thisC,dc); mdl_fun_h = str2func(gmdl);
        try
            mdl = fitnlm(xs',mSig',mdl_fun_h,x0,'options',statsetfitnlm);
        catch
            continue;
        end
        Rs(jj) = mdl.Rsquared.Ordinary;
%         aRs(jj) = mdl.Rsquared.Adjusted;
    end
    inds = Rs == max(Rs);
    if sum(inds) > 1
        inds = find(inds);
        coeffs = [pks(inds(1)),locs(inds(1)),w(inds(1))];
    end
    if sum(inds) == 1
        coeffs = [pks(inds),locs(inds),w(inds)];
    end
    if sum(inds) == 0
        coeffs = [max_mSig,loc_max_mSig,length(xs)/10];
    end
    [gmdl,x0] = buildGaussianModel(order,coeffs,dc); mdl_fun_h = str2func(gmdl);
    try
        mdl = fitnlm(xs',mSig',mdl_fun_h,x0,'options',statsetfitnlm);
    catch
        continue;
    end
    rasterF(ii,:) = mdl_fun_h(mdl.Coefficients.Estimate,xs);
    rs(ii,:) = [mdl.Coefficients.Estimate' mdl.Rsquared.Ordinary];
end
