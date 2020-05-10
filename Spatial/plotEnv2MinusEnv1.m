function plotEnv2MinusEnv1
ei = evalin('base','ei');
owsi = 1;
tei = ei(7:8);
allA_P = getRasters(tei,1:10);
% SI1 = find_spatial_information(tei,allA_P,owsi);
allAP_P = getRasters(tei,13:20);
% SI2 = find_spatial_information(tei,allAP_P,owsi);
dSIsM2 = findPearsonCorrelation(allA_P,allAP_P);
% dSIsM2 = SI2.SIm - SI1.SIm;
% dSIsM2 = SI2.M - SI1.M;


minVal = min([dSIsM2]);
maxVal = max([dSIsM2]);
bins = linspace(minVal,maxVal,20);
% bins = -30:5:70;
[bar1 xs] = hist(dSIsM2,bins);
% [bar2 xs] = hist(dSIsRSC,bins);
allBars = [100*bar1/sum(bar1)];%;100*bar2/sum(bar2)];

figure(100);clf;
hbars = bar(xs,allBars');
set(hbars(1),'facecolor','b');
% set(hbars(2),'facecolor','r');
xlim([bins(1) bins(end)]);
set(gca,'TickDir','out','FontSize',14,'FontWeight','Bold');
xlabel('Spatial Information Difference (bits)');
ylabel('Percentage');