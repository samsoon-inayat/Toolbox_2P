function raster = make_gauss_fit_raster(R,cn)

gfom = R.gauss_fit_on_mean;
xs = 1:size(R.sp_rasters1,2);
coeff = gfom.coefficients_Rs_trials(:,:,cn);
raster = NaN(size(R.sp_rasters1(:,:,1)));
for ii = 1:size(coeff,1)
    tcoeff = coeff(ii,:);
    fitplot = gauss_fit(xs,tcoeff(1:3),R.gauss_fit_on_mean.gauss1Formula);
    raster(ii,:) = fitplot;
end
% 
% plot(size(thisRaster,1)*fitplot/max(fitplot),'linewidth',0.5,'color','m');