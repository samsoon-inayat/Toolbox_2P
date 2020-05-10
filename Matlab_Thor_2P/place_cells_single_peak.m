function props = place_cells_single_peak (mSig,ccSignalA,varargin)

p = inputParser;
default_cm_per_bin = 157/100;
default_min_pc_width = 3;
default_max_pc_width = 100;

addRequired(p,'mSig',@isnumeric);
addRequired(p,'ccSignalA',@isnumeric);
addOptional(p,'cm_per_bin',default_cm_per_bin,@isnumeric);
addOptional(p,'max_pc_width',default_max_pc_width,@isnumeric);
addOptional(p,'min_pc_width',default_min_pc_width,@isnumeric);
parse(p,mSig,ccSignalA,varargin{:});
cm_per_bin = p.Results.cm_per_bin;
min_pc_width = p.Results.min_pc_width;
max_pc_width = p.Results.max_pc_width;

props.widths = NaN;
props.centers = NaN;
props.peaks = NaN;

xs = 1:length(mSig);
mSig = fixMSig(xs,mSig);
uxs = linspace(1,length(mSig),length(mSig)*10);
umSig = interp1(xs,mSig,uxs);
umSig = applyGaussFilt(umSig,10);

% skm = skewness(umSig);
% kurtm = kurtosis(umSig);
% 
% [p,tbl,stats] = kruskalwallis(ccSignalA,[],'off');
% if p>0.001
%     return;
% end
% threshold = mean(umSig);
% threshold = (0.2*(max(mSig) - min(mSig)));
% umSig(umSig<threshold) = 0;
bad = 1;
for ii = 1
    try
        cmdTxt = sprintf('[f,gof,output1] = fit(uxs'',umSig'',''gauss%d'');',ii); %[f1,gof1,output1] = fit(uxs',umSig','gauss1');
        eval(cmdTxt);
    catch
        break;
    end
    fy = f(uxs);
%     plotTheCurves(uxs,umSig,ccSignalA,fy,cm_per_bin);
    if gof.rsquare > 0.5
        bad = 0;
        break;
    else
        continue;
    end
end
if bad
    display('bad fitting');
    return;
end
fitOrder = ii;
titleTxt = sprintf('%.2f',gof.rsquare);
% plotTheCurves(uxs,umSig,ccSignalA,fy,cm_per_bin,titleTxt);
% if fitOrder > 1
%     combs = nchoosek(1:fitOrder,2);
% end

for jj = 1%:fitOrder
    cmdTxt = sprintf('as(jj) = f.a%d;bs(jj) = f.b%d;cs(jj) = f.c%d;',jj,jj,jj);
    eval(cmdTxt);
end

% find widths
PWs = 2.36*cs./sqrt(2)*cm_per_bin; % width at half max


% for ii = 1%:size(combs,1)
%     dbs(ii) = abs(cm_per_bin*(bs(combs(ii,2)) - bs(combs(ii,1))));
% end
% 
% if any(dbs<10)
%     n = 0;
% end

% apply rules
PWinds = PWs > min_pc_width & PWs < max_pc_width; % width has to be greater than 5 cm and less than 120 cm
bsInds = bs > 1 & bs < max(xs); % the center of place field has to be greater that 1 and less than the max of x (or belt)
asInds = as > (0.2*(max(mSig) - min(mSig))) & as > 10; % the amplitude threshold below which its noise
inds = PWinds & bsInds & asInds;
if ~any(inds)
    return;
end
props.widths = PWs(inds);
props.centers = bs(inds) * cm_per_bin;
props.peaks = as(inds);
% plotTheCurves(uxs,umSig,ccSignalA,fy,cm_per_bin);
n = 0;


function plotTheCurves(uxs,umSig,ccSignalA,fy,cm_per_bin,titleText)

[pks,locs,wpk,pkp] = findpeaks(umSig,'MinPeakDistance',50,'MinPeakProminence',50);
thresh = mean(umSig) + 0.5*std(umSig);
inds = find(pks > thresh);
pks = pks(inds);locs = locs(inds);wpk = wpk(inds);pkp = pkp(inds);
figure(10002);clf
subplot 221;imagesc(ccSignalA);set(gca,'XTick',(1:5:50),'XTickLabel',round(((1:5:50))*cm_per_bin));
subplot 222;imagesc(corrcoef(ccSignalA));
colorbar;
set(gca,'Ydir','Normal');
subplot(2,2,3);
plot(uxs*cm_per_bin,umSig,'k');hold on;
plot([uxs(1) uxs(end)]*cm_per_bin,[mean(umSig) mean(umSig)],'m');
if ~isnan(fy)
    plot(uxs*cm_per_bin,fy,'r');
end
plot(uxs(locs)*cm_per_bin,umSig(locs),'*r');
title(titleText);
set(gca,'XTick',round((uxs(1:50:500))*cm_per_bin));
xlim([uxs(1) uxs(end)]*cm_per_bin);
% title(sprintf('%.3f - %.3f',skm,kurtm));



%     umSig = interp1(xs,mSig,uxs);
%     umSig = applyGaussFilt(umSig,10);
    
%     [pks,locs,wpk,pkp] = findpeaks(umSig,'MinPeakDistance',50,'MinPeakProminence',50);
%     thresh = mean(umSig) + 0.5*std(umSig);
%     inds = find(pks > thresh);
%     pks = pks(inds);locs = locs(inds);wpk = wpk(inds);pkp = pkp(inds);
%     nPeaks(ii) = length(inds);
%     if plotFlag
%         figure(101);clf;
%         plot(xs,mSig);hold on;
%         plot(uxs,umSig);
%         plot(uxs(locs),pks,'*','color','r');
%         title(nPeaks(ii));
%         n = 0;
%     end
%     
%     mSig = applyGaussFilt(nanmean(raster),7);