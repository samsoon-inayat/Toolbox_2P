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
RsC = oC.Rs;% get_rasters_data(ei_C,selContexts,rasterNames);
RsA = oA.Rs;% get_rasters_data(ei_A,selContexts,rasterNames);
% typeP = {'all','vals'
%%
if 1

for rr = 1:size(RsC,1)
    numset1(rr) = length(RsC{rr,1}.iscell);
    numset2(rr) = length(RsA{rr,1}.iscell);
end
numset1(1) = (sum(RsC{1,1}.cns(:,2) == 1) + sum(RsC{1,1}.cns(:,2) == 2))/2;
[h,p,ci,stat] = ttest2(numset1,numset2)
[mC,semC] = findMeanAndStandardError(numset1);
[mA,semA] = findMeanAndStandardError(numset2);
    mVar = [mC mA];
    semVar = [semC semA];
    combs = [1 2]; 
    xdata = 1:2;
    tcolors
    hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.5 1],'color','w');
    hold on;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',5,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.5,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.01);
    set(gca,'xlim',[0.25 2.75],'ylim',[0 maxY],'FontSize',6,'FontWeight','Bold','TickDir','out','xcolor','k','ycolor','k');
    xticks = [1 2]; xticklabels = {'C-TG','A-TG'}; 
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
    for ii = 1:length(hbs)
        set(hbs(ii),'facecolor','none','edgecolor',tcolors{ii});
    end
%     xtickangle(30);
%     rectangle(gca,'Position',[0.75 30.5 1 3.5],'edgecolor','k','facecolor','k');
%     text(1.85,30.5,'Trials','FontSize',5);
%     rectangle(gca,'Position',[6 30.5 1 3.5],'edgecolor','k');
%     text(7.2,30.5,'Inter-Trials','FontSize',5);
    changePosition(gca,[0.1 0.13 -0.2 -0.1])
    put_axes_labels(gca,{[],[0 0 0]},{{'Number of Cells'},[0 0 0]});
    
    save_pdf(hf,mData.pdf_folder,'Figure_Number_Of_cells.pdf',600);
return;
end

