function figure_distribution_of_zMI

mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size;
ei_C = evalin('base','ei10_C2'); 
ei_A = evalin('base','ei10_A2'); 


selContexts = [1 2 3 4];
rasterNames = {'airD','airD','airD','airD'};

RsC = get_rasters_data(ei_C,selContexts,rasterNames);
mRsC = calc_mean_rasters(RsC,1:10);
RsC = find_responsive_rasters(RsC,1:10);
% view_population_vector(Rs,mRs,300);
[resp_fractionC,resp_valsC,OIC,mean_OIC] = get_responsive_fraction(RsC)

RsA = get_rasters_data(ei_A,selContexts,rasterNames);
mRsA = calc_mean_rasters(RsA,1:10);
RsA = find_responsive_rasters(RsA,1:10);
% view_population_vector(Rs,mRs,400);
[resp_fractionA,resp_valsA,OIA,mean_OIA] = get_responsive_fraction(RsA)

n = 0;
%%

tcolors = {'k','r'};
distD(:,1) = out_C.allVals_an';
distD(:,2) = out_A.allVals_an';
[distDo,allVals,allValsG] = plotDistributions(distD);
minBin = min(allVals);
maxBin = max(allVals);
[h,p,ks2stat] = kstest2(allValsG{1},allValsG{2});
[h,p,cd,ks2stat] = ttest2(allValsG{1},allValsG{2});
%%
incr = 0.001; %maxBin =
hf = figure(8);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.85 1],'color','w');
hold on;
%    [ha,hb,hca] = plotDistributions(distD,'colors',tcolors,'maxY',maxBin,'cumPos',[0.5 0.26 0.25 0.5],'min',minBin,'incr',incr,'max',maxBin);
[ha,hb,hca] = plotDistributions(distDo,'colors',tcolors,'maxY',100,'min',minBin,'incr',incr,'max',maxBin);
set(gca,'FontSize',6,'FontWeight','Bold','TickDir','out','xcolor','k','ycolor','k');
changePosition(gca,[0.09 0.13 -0.05 -0.13]);
put_axes_labels(gca,{'Average Firing Rate (Hz)',[0 0 0]},{{'Percentage','of Neurons'},[0 0 0]});
save_pdf(hf,mData.pdf_folder,sprintf('Distribution_firing_rate'),600);