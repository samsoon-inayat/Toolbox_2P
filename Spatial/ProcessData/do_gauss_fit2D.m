function [rasterF,mdl,coeff_rs] = do_gauss_fit2D(raster,statsetfitnlm,orderdc)
order = orderdc(1);
dc = orderdc(2);
rasterF = NaN(size(raster));
coeff_rs = NaN(size(raster,1),1+order*5);

xs = 1:size(raster,2);
ys = 1:size(raster,1);
[x,y] = meshgrid(xs,ys);
mval = nanmean(raster(:));
% tbl = table(x(:),y(:),raster(:));
pred = [x(:) y(:)];
mSig = raster(:);
[max_mSig,loc_max_mSig] = max(mSig);
[pks,locs,w,p] = findpeaks(mSig);
inds = p > nanmean(p);
if sum(inds) > 1
    pks = pks(inds);locs = locs(inds);w = w(inds);
end
Rs = NaN(1,length(pks));
for jj = 1:length(pks)
    locx = pred(locs(jj),1); locy = pred(locs(jj),2);
    thisC = [pks(jj),locx,w(jj),locy,w(jj)];
    [gmdl,x0] = buildGaussianModel2D(order,thisC,dc); mdl_fun_h = str2func(gmdl);
    if dc
        x0(end) = mval;
    end
    try
        mdl = fitnlm(pred,mSig,mdl_fun_h,x0,'options',statsetfitnlm);
    catch
        continue;
    end
    Rs(jj) = mdl.Rsquared.Ordinary;
%         aRs(jj) = mdl.Rsquared.Adjusted;
end
inds = Rs == max(Rs);
if sum(inds) > 1
    inds = find(inds);
    jj = inds(1);
    locx = pred(locs(jj),1); locy = pred(locs(jj),2);
    coeffs = [pks(jj),locx,w(jj),locy,w(jj)];
end
if sum(inds) == 1
    jj = inds;
    locx = pred(locs(jj),1); locy = pred(locs(jj),2);
    coeffs = [pks(jj),locx,w(jj),locy,w(jj)];
end
if sum(inds) == 0
    locx = pred(loc_max_mSig,1); locy = pred(loc_max_mSig,2);
    coeffs = [max_mSig,locx,length(xs)/10,locy,length(xs)/10];
end
[gmdl,x0] = buildGaussianModel2D(order,coeffs,dc); 
mdl_fun_h = str2func(gmdl);
mdl = fitnlm(pred,raster(:),mdl_fun_h,x0,'options',statsetfitnlm);

rasterF = mdl_fun_h(mdl.Coefficients.Estimate,pred);
rasterF = reshape(rasterF,size(raster,1),size(raster,2));
coeff_rs = [mdl.Coefficients.Estimate' mdl.Rsquared.Ordinary];

