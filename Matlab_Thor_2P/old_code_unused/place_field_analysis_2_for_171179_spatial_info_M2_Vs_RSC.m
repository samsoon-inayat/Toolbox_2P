function place_field_analysis_2_for_171179_spatial_info_M2_Vs_RSC
ei = evalin('base','ei');
owsi = 0;
tei = ei(1:2);
allA_P = getRasters(tei,1:20);
% SI1 = find_spatial_information(tei,allA_P,owsi);
MI1 = find_mutual_information(tei,allA_P,owsi);
allAP_P = getRasters(tei,13:23);
MI2 = find_mutual_information(tei,allAP_P,owsi);
% SI2 = find_spatial_information(tei,allAP_P,owsi);
% dSIsM2 = findPearsonCorrelation(allA_P,allAP_P);
% dSIsM2 = MI2.zMI;
dSIsM2 = MI2.zMI - MI1.zMI;
dSIsM2 = MI2.zMI;
% dSIsM2 = SI2.M - SI1.M;

tei = ei(3:4);

allA_P = getRasters(tei,1:11);
MI1 = find_mutual_information(tei,allA_P,owsi);
allAP_P = getRasters(tei,13:23);
MI2 = find_mutual_information(tei,allAP_P,owsi);
% SI2 = find_spatial_information(tei,allAP_P,owsi);
% dSIsRSC = findPearsonCorrelation(allA_P,allAP_P);
% dSIsRSC = MI2.zMI;
dSIsRSC = MI2.zMI - MI1.zMI;
dSIsRSC = MI2.zMI;
% dSIsRSC = SI2.SIm - SI1.SIm;
% dSIsRSC = SI2.M - SI1.M;


minVal = min([dSIsM2 dSIsRSC]);
maxVal = max([dSIsM2 dSIsRSC]);
% maxVal = 2;
bins = linspace(minVal,maxVal,20);
% bins = -30:5:70;
[bar1 xs] = hist(dSIsM2,bins);
[bar2 xs] = hist(dSIsRSC,bins);
allBars = [100*bar1/sum(bar1);100*bar2/sum(bar2)];

figure(100);clf;
hbars = bar(xs,allBars');
set(hbars(1),'facecolor','b');
set(hbars(2),'facecolor','r');
xlim([bins(1) bins(end)]);
set(gca,'TickDir','out','FontSize',14,'FontWeight','Bold');
xlabel('Mutual Information Difference (bits)');
ylabel('Percentage');
% save2pdf('SI_diff_M2_RSC.pdf',gcf,600);
[p,h] = ranksum(dSIsM2,dSIsRSC,'Tail','left')


% plotMarkers(ei1,ei1.b.air_puff_r(1:23),[],101);
% xlim([550 1000]); ylim([-0.2 1]);

n = 0;