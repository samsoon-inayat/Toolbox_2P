function out =  getFits(thispFolder,contextNumber,stimMarker,Rs,thisRasterType,trials,ow)

fileName = makeName(sprintf('gauss_fit_on_means_%d_%s_%s.mat',contextNumber,stimMarker,thisRasterType),thispFolder);
if exist(fileName,'file') && ow == 0
    out = load(fileName);
    return;
end
hWaitBar = waitbar(0,sprintf('Finding gaussian fits ... plz wait'));


rasters = Rs.sp_rasters_nan_corrected(trials,:,:);


numCells = size(rasters,3);
xs = 1:size(rasters,2);

uxs = linspace(1,length(xs),length(xs)*10);

gauss1Formula = '';
rs = NaN(numCells,1);
adj_rs = rs;
worked = rs;
coeffs = NaN(numCells,3);
try

umSig = NaN(numCells,length(uxs));
parfor ii = 1:numCells
    mSig = nanmean(rasters(:,:,ii));
    umSig(ii,:) = interp1(xs,mSig,uxs);
    umSig(ii,:) = applyGaussFilt(umSig(ii,:),10);
end

[f,~,~] = fit(uxs',umSig(1,:)','gauss1');
gauss1Formula = formula(f);

tic
parfor ii = 1:numCells
%     waitbar(ii/numCells,hWaitBar,sprintf('Gaussian Fitting - Processing cell %d/%d',ii,numCells));
%     raster = rasters(:,:,ii);
%     mSig = nanmean(raster);
    try
%         mSig = fixMSig(xs,mSig);
%         umSig = interp1(xs,mSig,uxs);
%         umSig = applyGaussFilt(umSig,10);
        [f,gof,~] = fit(uxs',umSig(ii,:)','gauss1');
    catch
        worked(ii) = 0;
%         disp('Gaussian fitting error');
        continue;
    end
    worked(ii) = 1;
%     allgof(ii) = gof;
%     alloutput(ii) = output1;
    rs(ii) = gof.rsquare;
    adj_rs(ii) = gof.adjrsquare;
    coeffs(ii,:) = coeffvalues(f);
end
catch
    disp(sprintf('Error occurred in Context %d --- %s -- %s',contextNumber,stimMarker,thisRasterType));
end
close(hWaitBar);
toc
out.gauss1Formula = gauss1Formula;
out.worked = worked;
out.rs = rs;
out.adj_rs = adj_rs;
out.coeffs = coeffs;
save(fileName,'-struct','out','-v7.3');
