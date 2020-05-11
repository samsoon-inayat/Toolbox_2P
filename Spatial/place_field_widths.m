function place_field_widths


% place_cell_properties.m




ei = evalin('base','ei');

zMITh = 1:11;
recs = {[1 2],[3 4],[5 6],[7 8]};
for zz = 1:length(zMITh)
    for ii = 1:4
        normSel = 1;
        rec = recs{ii};
        pcw = [];ind = 1;
        allRVals = [];
        for jj = 1:2
            ptc = getPositionTuning(ei{rec(jj)},'sp');
            zMIs = ei{rec(jj)}.rasters.zMI;
            cois = find(zMIs>zMITh(zz));
            if isempty(cois)
                continue;
            end
            ptc = ptc(cois,:);
%             allRVals = [allRVals ei{rec(jj)}.rasters.rasters(:)'];
            for kk = 1:size(ptc,1)
                mSig = ptc(kk,:);
                raster = ei{rec(jj)}.rasters.rasters(:,:,cois(kk));
                if max(raster(:))<100%max(ei{rec(jj)}.deconv.caSigAll{kk})<100%  
                    continue;
                end
                pcw(ind) = place_cell_properties(mSig,raster,'cm_per_bin',142/50);
                ind = ind + 1;
            end
            n=0;
        end
        allpcw(zz,ii) = nanmean(pcw);
%         pcwTemp = pcw(~isnan(pcw));
%         allpcwsem(zz,ii) = std(pcwTemp)/sqrt(length(pcwTemp));
        allpcwVals{zz,ii} = pcw;
    end
end
nCells(:,1) = nanmean(allpcw(:,[1 2]),2);
nCells(:,2) = nanmean(allpcw(:,[3 4]),2);

figure(1001);clf;%tempPos = get(gcf,'Position');tempPos(3:4) = [4.5 4.5];set(gcf,'Units','Inches');set(gcf,'Position',tempPos);
hp = plot(nCells,'linewidth',2);
cols = {'b','r'};
for ii = 1:2
    set(hp(ii),'Color',cols{ii});
end
xlabel('Normalized Mutual Information Threshold (Z-Score)');
ylabel('Mean Place Cell Width (cm)');
set(gca,'linewidth',1.5,'FontSize',16,'FontWeight','Bold');
box off;

hold on;
legendText = {'APP(NL-G-F)-Thy1-GCaMP','Thy1-GCaMP'};
thisCols =  {'r','b'};
x1 = 1; x2 = x1+1; y1 = (16:1:100); y1 = y1(1:4); y2 = y1;
legendFontSize = 16;
legs = legendText;
for ii = 1:length(legs)
    plot([x1 x2],[y1(ii) y2(ii)],'color',thisCols{ii},'linewidth',2);
    text(x2+0.5,y1(ii),sprintf('%s',legs{ii}),'Color',thisCols{ii},'FontSize',legendFontSize);
end
save2pdf('zscoreAnalysis.pdf',gcf,600);



%%
selZ = 5;
pcws_c = [allpcwVals{selZ,1} allpcwVals{selZ,2}];
pcws_d = [allpcwVals{selZ,3} allpcwVals{selZ,4}];

pcws_c(pcws_c>120) = [];pcws_c(pcws_c<5) = [];
pcws_d(pcws_d>120) = [];pcws_d(pcws_d<5) = [];

pcws_c(pcws_c==NaN) = [];
pcws_d(pcws_d==NaN) = [];

vals1 = pcws_c;
vals2 = pcws_d;

minVal = min([vals1 vals2]);
maxVal = max([vals1 vals2]);
% maxVal = 2;
bins = linspace(minVal,maxVal,20);
% bins = -30:5:70;
[bar1 xs] = hist(vals1,bins);
[bar2 xs] = hist(vals2,bins);
allBars = [100*bar1/sum(bar1);100*bar2/sum(bar2)];

figure(10001);clf;
hbars = bar(xs,allBars');
set(hbars(1),'facecolor','b');
set(hbars(2),'facecolor','r');
xlim([bins(1) bins(end)]);
set(gca,'TickDir','out','FontSize',14,'FontWeight','Bold','linewidth',1.5);box off;
xlabel('Place Field Widths (cm)');
ylabel('Percentage');
hold on;
legendText = {'APP(NL-G-F)-Thy1-GCaMP','Thy1-GCaMP'};
legsN = [2 1];
thisCols = {'r','b'};
x1 = 15; x2 = x1+5; y1 = (37:2:100); y1 = y1(1:2); y2 = y1;
legendFontSize = 11;
legs = legendText;
for ii = 1:length(legs)
    plot([x1 x2],[y1(ii) y2(ii)],'color',thisCols{ii},'linewidth',2);
    eval(sprintf('len = length(vals%d);',legsN(ii)));
    text(x2+0.5,y1(ii),sprintf('%s, (n=%d cells, 2 animals)',legs{ii},len),'Color',thisCols{ii},'FontSize',legendFontSize);
end
xlim([0 max(xs)]);
text(5,44,sprintf('For cells with Normalized Mutual Information (Z-Score) > 5'),'FontSize',14,'FontWeight','bold');

axesPosO = get(gca,'Position'); 
axesPos = get(gca,'Position');
cAllBars = cumsum(allBars,2);
subAxesPos(1) = axesPos(1) + 0.27;
subAxesPos(2) = axesPos(2) + 0.25;
subAxesPos(3) = 0.15;
subAxesPos(4) = 0.25;
axes('Position',subAxesPos);
plot(xs,cAllBars(1,:),'linewidth',1.5,'color','b');hold on;
plot(xs,cAllBars(2,:),'linewidth',1.5,'color','r');
set(gca,'TickDir','out','FontSize',12,'FontWeight','Bold','linewidth',1.5);box off
ht = title('Cumulative Distribution');
set(ht,'FontSize',12);
tPos = get(ht,'Position');
set(ht,'Position',(tPos + [-10 3 0]));
xlim([0 max(xs)]);

subAxesPos(1) = axesPosO(1) + 0.54;
subAxesPos(2) = axesPosO(2) + 0.25;
subAxesPos(3) = 0.2;
subAxesPos(4) = 0.25;
axes('Position',subAxesPos);
hp = plot(nCells,'linewidth',2);
cols = {'b','r'};
for ii = 1:2
    set(hp(ii),'Color',cols{ii});
end
xlabel('Norm. MI Thresh(Z-Score)');
ylabel('Mean PF Width (cm)');
set(gca,'TickDir','out','FontSize',12,'FontWeight','Bold','linewidth',1.5);box off
xlim([0 11]);


% save2pdf('SI_diff_M2_RSC.pdf',gcf,600);
[p,h] = ranksum(vals1,vals2,'Tail','left')
[ht,pt] = ttest2(vals1,vals2,'Tail','left')
[cM cE cS] = findMeanAndStandardError(pcws_c);
[dM dE dS] = findMeanAndStandardError(pcws_d);

% text(0,-40,sprintf('***, p = %.3f, Student''s t-test',p),'FontSize',12,'FontWeight','bold');

figure(10001);
save2pdf('place_field_widths.pdf',gcf,600);
n = 0;

%%

