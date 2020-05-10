function [mdl,result] = fit2DGauss(thisRaster)

[ydim,xdim] = size(thisRaster);

xx = 1:xdim;
yy = 1:ydim;

[X,Y] = meshgrid(xx,yy);

xData = reshape(X,numel(X),1);
yData = reshape(Y,numel(Y),1);
zData = reshape(thisRaster,numel(X),1);
xyData = [xData,yData];

statsetfitnlm = statset('fitnlm');
statsetfitnlm.MaxIter = 1000;
statsetfitnlm.TolFun = 1e-8;
% statsetfitnlm.Display = 'iter';
statsetfitnlm.TolX = 1e-8;
statsetfitnlm.UseParallel = 1;

coeffs = [1,floor(xdim)/2,floor(xdim)/10,floor(ydim)/2,floor(ydim)/10];

mdl = fitnlm(xyData,zData,@modelfun_gaussian2D,[0,coeffs,coeffs,coeffs],'options',statsetfitnlm);

result = modelfun_gaussian2D(mdl.Coefficients.Estimate,xyData);
result = reshape(result,ydim,xdim);
                                                        
% [coeffs,resnorm,residual,exitflag,output,lambda,jacobian] = lsqcurvefit(@modelfun_gaussian2D,[0,1,floor(xdim)/2,floor(xdim)/10,floor(ydim)/2,...
%                                                             floor(ydim)/10],xyData,zData,[],[],options);
% 
% rF = modelfun_gaussian2D(coeffs,xyData);
% rF = reshape(rF,ydim,xdim);
% result = rF;
% 
% figure(1000);clf;
% subplot 121; imagesc(thisRaster);
% subplot 122; imagesc(rF);
                                                        
n = 0;