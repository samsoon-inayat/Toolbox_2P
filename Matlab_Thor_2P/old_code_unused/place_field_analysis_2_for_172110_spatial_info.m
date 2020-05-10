function place_field_analysis_2_for_172110_spatial_info
ei{1} = evalin('base','ei{5}');
ei{2} = evalin('base','ei{6}');

% plotMarkers(ei1,ei1.b.air_puff_r(1:23),[],101);

allA_P = getRasters(ei,1:11,1);
allAP_P = getRasters(ei,13:23,1);


bins = -0.5:0.25:2.5;
[bar1 xs] = hist(dSIsM2,bins);
[bar2 xs] = hist(dSIsRSC,bins);
allBars = [100*bar1/sum(bar1);100*bar2/sum(bar2)];

figure(100);clf;bar(xs,allBars');
xlim([bins(1) bins(end)]);
set(gca,'TickDir','out','FontSize',14,'FontWeight','Bold');
xlabel('Spatial Information Difference (bits)');
ylabel('Percentage');
save2pdf('SI_diff_M2_RSC.pdf',gcf,600);



plotMarkers(ei2,ei1.b.air_puff_r(1:23),[],101);
% xlim([550 1000]); ylim([-0.2 1]);

n = 0;