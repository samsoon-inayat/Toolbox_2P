function light_figure

mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size;
ei = evalin('base','ei'); 

selContexts = [1 4 6];
rasterNames = {'light22T','light22T','light22T'};
Rs = get_rasters_data(ei,selContexts,rasterNames);
Rs = find_responsive_rasters(Rs,1:10);
[resp_fractionL,resp_valsL,OIL,mean_OIL,resp_ORL,resp_OR_fractionL,resp_ANDL,resp_AND_fractionL] = get_responsive_fraction(Rs);

selContexts = [2 7];
rasterNames = {'air55T','air55T'};
Rs = get_rasters_data(ei,selContexts,rasterNames);
Rs = find_responsive_rasters(Rs,1:10);
[resp_fractionA,resp_valsA,OIA,mean_OIA,resp_ORA,resp_OR_fractionA,resp_ANDA,resp_AND_fractionA] = get_responsive_fraction(Rs);


selContexts = [3 4 5];
rasterNames = {'airD','airD','airD'};
% Rs = get_rasters_data(ei_11_15,selContexts,rasterNames);
Rs = get_rasters_data(ei,selContexts,rasterNames);
Rs = find_responsive_rasters(Rs,1:10);
mR = calc_mean_rasters(Rs,1:10);
[CRc,aCRc,mRR] = find_population_vector_corr(Rs,mR,1);
[resp_fractionC,resp_valsC,OIC,mean_OIC,resp_ORC,resp_OR_fractionC,resp_ANDC,resp_AND_fractionC] = get_responsive_fraction(Rs);
for ii = 1:length(resp_valsC)
    selResp{ii} = logical(resp_valsC{ii}(:,1));
end
[all_corr_C,all_corr_cell_C,mean_corr_C,mean_cell_corr_C,xs_C] = find_population_vector_corr_remap(Rs,mR,selResp);

n = 0;

%%
if 0
% out_CC = all_out_C.all_corr_an{3};
% out_AA = all_out_A.all_corr_an{1};
    temp = mean_cell_corr_C(:,:,1);
    mask = ones(size(temp)); mask = triu(mask,1) & ~triu(mask,2);
    for ii = 1:size(mean_cell_corr_C,3)
        temp = mean_cell_corr_C(:,:,ii);
        var_C(ii,:) = temp(mask)';
    end
    dataT = array2table([var_C]);
    dataT.Properties.VariableNames = {'C12','C23'};
    colVar1 = [1 2];    
    within = table(colVar1');
    within.Properties.VariableNames = {'Condition'};
    within.Condition = categorical(within.Condition);
    ra = repeatedMeasuresAnova(dataT,within);

    mVar = ra.est_marginal_means.Mean;semVar = ra.est_marginal_means.Formula_StdErr;
    combs = ra.mcs.combs; p = ra.mcs.p; h = ra.mcs.p < 0.05;
    xdata = [1 2];
    colors = mData.colors;
    hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.5 1],'color','w');
    hold on;
    tcolors ={colors{1};colors{2};colors{3};colors{1};colors{2};colors{3};};
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.2,'sigTestName','','sigLineWidth',0.25,'BaseValue',0,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.7,'sigLinesStartYFactor',0.1);
    % plot([0.5 11],[-0.5 0.5],'linewidth',1.5)
    set(gca,'xlim',[0.25 xdata(end)+0.75],'ylim',[0 maxY],'FontSize',6,'FontWeight','Bold','TickDir','out');
    xticks = xdata(1:end)+0; xticklabels = {'C3-C4','C4-C3'''};
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
    xtickangle(30);
    changePosition(gca,[0.06 0.03 0.02 -0.11]);
    put_axes_labels(gca,{[],[0 0 0]},{{'Correlation'},[0 0 0]});

    save_pdf(hf,mData.pdf_folder,'cell corr remap',600);
return;
end

%%

    sel_out = mean_corr_C;
    xs = xs_C;


ff = makeFigureRowsCols(107,[1 0.5 2 2],'RowsCols',[3 3],...
    'spaceRowsCols',[0.05 0.07],'rightUpShifts',[0.071 0.1],'widthHeightAdjustment',...
    [-95 -95]);
set(gcf,'color','w');
set(gcf,'Position',[1 6 3.4 3.4]);
FS = mData.axes_font_size;
maskdisp = triu(ones(4,4),0);

for rr = 1:3
    for cc = 1:3
        if ~maskdisp(rr,cc)
            delete(ff.h_axes(rr,cc));
            continue;
        end
        axes(ff.h_axes(rr,cc));
        imagesc(sel_out{rr,cc});
        box off;
        set(gca,'Ydir','Normal','linewidth',1,'FontSize',FS-1,'FontWeight','Bold');
        if rr == cc
            cols = size(sel_out{rr,cc},2);
            colsHalf = round(cols/2);
            ts = round(xs);
            set(gca,'YTick',[1 colsHalf cols],'YTickLabel',[0 ts(colsHalf) ts(cols)]);
            h = xlabel('Position (cm)');%    changePosition(h,[0 0 0]);
            h = ylabel('Position (cm)');%    changePosition(h,[0 0 0]);
            cols = size(sel_out{rr,cc},2);
            colsHalf = round(cols/2);
            ts = round(xs);
            set(gca,'XTick',[1 colsHalf cols],'XTickLabel',[0 ts(colsHalf) ts(cols)]);

        else
            axis off
        end
        minC = min(sel_out{rr,cc}(:));
        maxC = max(sel_out{rr,cc}(:));
        hc = putColorBar(ff.h_axes(rr,cc),[0.0 0.03 0 -0.05],[minC maxC],5,'eastoutside',[0.07 0.08 0.1 0.13]);
    end
end
save_pdf(ff.hf,mData.pdf_folder,sprintf('pos corr remap'),600);
return;