function out =  getFits_myGaussFit_1(thispFolder,contextNumber,stimMarker,Rs,thisRasterType,trials,ow)

fileName = makeName(sprintf('gauss_fit_on_means_1_%d_%s_%s.mat',contextNumber,stimMarker,thisRasterType),thispFolder);
if exist(fileName,'file') && ow == 0
    out = load(fileName);
    return;
end
% hWaitBar = waitbar(0,sprintf('Finding gaussian fits ... plz wait'));

rasters = Rs.sp_rasters_nan_corrected;%(trials,:,:);

numCells = size(rasters,3);
xs = 1:size(rasters,2);

% uxs = linspace(1,length(xs),length(xs)*100);

worked = zeros(numCells,1);
% rs = NaN(numCells,1);
% mSigA = NaN(numCells,length(xs));
% umSig = NaN(numCells,length(uxs));
% umSigF = NaN(size(umSig));

% parfor ii = 1:numCells
%     mSig = nanmean(rasters(:,:,ii));
%     mSigA(ii,:) = applyGaussFilt(mSig,4);
% %     temp = interp1(xs,mSig,uxs);
% %     umSig(ii,:) = applyGaussFilt(temp,300);
% %     figure(1000);clf;plot(xs,mSig(ii,:),'b');hold on;plot(uxs,umSig(ii,:),'r');plot(uxs,temp,'g');
% end

statsetfitnlm = statset('fitnlm');
statsetfitnlm.MaxIter = 1000;
statsetfitnlm.TolFun = 1e-10;
% statsetfitnlm.Display = 'iter';
statsetfitnlm.TolX = statsetfitnlm.TolFun;
statsetfitnlm.UseParallel = 1;
% statsetfitnlm.RobustWgtFun = 'welsch';

% coefficients = NaN(4,4,numCells);
% orders = 2;%:-1:3;
coeffsrsM = NaN(numCells,4);
coeffsrsR = NaN(size(rasters,1),4,numCells);
tic
pctRunOnAll warning off;
parfor ii = 1:numCells
% for ii = 1:numCells
    try
        thisRaster = rasters(trials,:,ii);
        mSigO = nanmean(thisRaster);
        [~,~,coeffsrsM(ii,:)] = do_gauss_fit(xs,mSigO,statsetfitnlm,[1 0]);
        [~,~,coeffsrsR(:,:,ii)] = do_gauss_fit(xs,rasters(:,:,ii),statsetfitnlm,[1 0]);
%         rasterF = thisRaster;
%         for jj = 1:length(orders)
%             [rasterF,~,rsjj{jj}] = do_gauss_fit(xs,rasterF,statsetfitnlm,[orders(jj) 0]);
%         end
% %         weights = repmat(rsjj{jj}',1,size(rasterF,2));
%         mSig = nanmean(rasterF);
%         [mSigF,mdl] = do_gauss_fit(xs,mSig,statsetfitnlm,[1 0]);
%         [f,gof,outF] = fit(xs',mSigO','gauss1');
%         mSigF = feval(f,xs);
%         rst = [mdlO1.Rsquared.Ordinary mdl.Rsquared.Ordinary];
% %         rst = [mdlO1.Rsquared.Ordinary mdlO2.Rsquared.Ordinary mdlO3.Rsquared.Ordinary];
%         figure(1000);clf;
%         subplot 121;imagesc(thisRaster);colorbar;title(ii);set(gca,'YDir','Normal');
%         subplot 122;plot(rasterF');colorbar;title(ii);set(gca,'YDir','Normal');
%         subplot 122;
%         plot(xs,mSigO,'b');hold on;
%         plot(xs,mSigOF,'r','linewidth',1.5);
%         plot(xs,mSigF,'g','linewidth',1.5);
% %         plot(xs,mSig,'m','linewidth',1.5);
%         title(rst);
%         rsjj{1}
    catch
        continue;
    end
    worked(ii) = 1;
end
pctRunOnAll warning on;

% close(hWaitBar);
toc
out.gauss1Formula = [1 0];
out.worked = logical(worked);
out.coefficients_Rs_mean = coeffsrsM;
out.coefficients_Rs_trials = coeffsrsR;
save(fileName,'-struct','out','-v7.3');


%         [mdl2D,result] = fit2DGauss(thisRaster);
% %         coefficients(:,:,ii) = mdl.Coefficients.Variables;
%         rs(ii,:) = [mdl.Rsquared.Ordinary];