function figure1_number_PCs_threshold_comparison

mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size;
ei_C = evalin('base','ei10_C1'); 
ei_A = evalin('base','ei10_A1'); 
tcolors = colors;
tcolors = {'k','r'};

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

%%
for rr = 1:size(RsC,1)
    naC = logical(ones(length(RsC{rr,1}.info_metrics.ShannonMI_Zsh),1));
    naA = logical(ones(length(RsA{rr,1}.info_metrics.ShannonMI_Zsh),1));
    for cc = 1:size(RsC,2)
%         numset1(rr,cc) = 100*sum(isnan(RsC{rr,cc}.info_metrics.ShannonMI_Zsh))/length(RsC{rr,cc}.iscell);
%         numset2(rr,cc) = 100*sum(isnan(RsA{rr,cc}.info_metrics.ShannonMI_Zsh))/length(RsC{rr,cc}.iscell);
        numset1(rr,cc) = 100*sum(isnan(RsC{rr,cc}.info_metrics.ShannonMI_Zsh))/sum(RsC{rr,cc}.iscell);
        numset2(rr,cc) = 100*sum(isnan(RsA{rr,cc}.info_metrics.ShannonMI_Zsh))/sum(RsA{rr,cc}.iscell);
        naC = naC & isnan(RsC{rr,cc}.info_metrics.ShannonMI_Zsh');
        naA = naA & isnan(RsA{rr,cc}.info_metrics.ShannonMI_Zsh');
    end
    all_naC(rr,1) = 100*sum(naC)/length(naC);
    all_naA(rr,1) = 100*sum(naA)/length(naA);
end


%%
dataT = array2table([[1;1;1;1;1;2;2;2;2;2] [numset1;numset2]]);
dataT.Properties.VariableNames = {'Group','C1','C2','C3','C4'};
within = array2table([1 2 3 4]');
within.Properties.VariableNames = {'Cond'};
within.Cond = categorical(within.Cond);
dataT.Group = categorical(dataT.Group);
ra = repeatedMeasuresAnova(dataT,within,0.05);
% writetable(dataT,fullfile(mData.pdf_folder,'Figure_Number_Of_cells_zMI_NaN.xls'));
%%
mVar = ra.est_marginal_means.Mean;
semVar = ra.est_marginal_means.Formula_StdErr;
combs = ra.mcs.combs; p = ra.mcs.p; h = 2*ra.mcs.p < 0.05;
xdata = [1 2 3 4 6 7 8 9]; 

hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.9 1],'color','w');
hold on;
tcolors = {colors{1};colors{2};colors{3};colors{4};colors{1};colors{2};colors{3};colors{4}};
% tcolors = repmat(tcolors,2,1)';
[hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
    'ySpacing',5,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
    'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.01);

set(gca,'xlim',[0.25 9.75],'ylim',[0 maxY],'FontSize',6,'FontWeight','Bold','TickDir','out','xcolor','k','ycolor','k');
xticks = xdata; xticklabels = {'C1','C2','C3','C4'}; xticklabels = repmat(xticklabels,1,2);
set(gca,'xtick',xticks,'xticklabels',xticklabels);

for ii = 5:length(hbs)
    set(hbs(ii),'facecolor','none','edgecolor',tcolors{ii});
end
%     rectangle(gca,'Position',[0.75 30.5 1 3.5],'edgecolor','k','facecolor','k');
%     text(1.85,30.5,'Trials','FontSize',5);
%     rectangle(gca,'Position',[6 30.5 1 3.5],'edgecolor','k');
%     text(7.2,30.5,'Inter-Trials','FontSize',5);
changePosition(gca,[0.07 0.02 -0.01 -0.011])
put_axes_labels(gca,{[],[0 0 0]},{{'Percentage of','Silent Neurons'},[0 0 0]});

save_pdf(hf,mData.pdf_folder,'Figure_Number_Of_cells_zMI_NaN.pdf',600);
