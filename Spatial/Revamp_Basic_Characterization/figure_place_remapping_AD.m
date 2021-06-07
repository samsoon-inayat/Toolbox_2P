function figure_place_remapping_AD(fn,allRs,ccs)

mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size;
ei_C = evalin('base','ei10_C2'); 
ei_A = evalin('base','ei10_A2'); 

selContexts = [1 2 3 4];
rasterNames = {'airD','airD','airD','airD'};

RsC = get_rasters_data(ei_C,selContexts,rasterNames);
mRsC = calc_mean_rasters(RsC,1:10);
RsC = find_responsive_rasters(RsC,1:10);
[CR_C,aCR_C] = find_population_vector_corr(RsC,mRsC,100);
[resp_fractionC,resp_valsC,OIC,mean_OIC,resp_ORC,resp_OR_fractionC,resp_ANDC,resp_AND_fractionC] = get_responsive_fraction(RsC);
[all_corr_C,all_corr_cell_C,mean_corr_C,mean_cell_corr_C,xs_C] = find_population_vector_corr_remap(RsC,mRsC,resp_ORC);


n = 0;
%%
if 1
% out_CC = all_out_C.all_corr_an{3};
% out_AA = all_out_A.all_corr_an{1};
    temp = mean_cell_corr_C(:,:,1);
    mask = ones(size(temp)); mask = triu(mask,1) & ~triu(mask,2);
    for ii = 1:size(mean_cell_corr_C,3)
        temp = mean_cell_corr_C(:,:,ii);
        var_C(ii,:) = temp(mask)';
    end
    for ii = 1:size(mean_cell_corr_A,3)
        temp = mean_cell_corr_A(:,:,ii);
        var_A(ii,:) = temp(mask)';
    end
    dataT = array2table([[ones(size(var_C,1),1);2*ones(size(var_A,1),1)] [var_C;var_A]]);
    dataT.Properties.VariableNames = {'Group','C12','C23','C34'};
    dataT.Group = categorical(dataT.Group);
    colVar1 = [1 2 3];    
    within = table(colVar1');
    within.Properties.VariableNames = {'Condition'};
    within.Condition = categorical(within.Condition);
    ra = repeatedMeasuresAnova(dataT,within);
    
    dataTC = dataT(1:3,2:end);
    raC = repeatedMeasuresAnova(dataTC,within);
    
    dataTA = dataT(4:end,2:end);
    raA = repeatedMeasuresAnova(dataTA,within);

    mVar = ra.est_marginal_means.Mean;semVar = ra.est_marginal_means.Formula_StdErr;
    combs = ra.mcs.combs; p = ra.mcs.p; h = ra.mcs.p < 0.05;
    xdata = [1 2 3 5:7];
    colors = mData.colors;
    hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.5 1],'color','w');
    hold on;
    tcolors ={colors{1};colors{2};colors{3};colors{1};colors{2};colors{3};};
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.2,'sigTestName','','sigLineWidth',0.25,'BaseValue',0,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.7,'sigLinesStartYFactor',0.1);
    for ii = 4:length(hbs)
        set(hbs(ii),'facecolor','none','edgecolor',tcolors{ii});
    end
    % plot([0.5 11],[-0.5 0.5],'linewidth',1.5)
    set(gca,'xlim',[0.25 xdata(end)+0.75],'ylim',[0 maxY],'FontSize',6,'FontWeight','Bold','TickDir','out');
    xticks = xdata(1:end)+0; xticklabels = {'C12','C23','C34','C12','C23','C34'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
    xtickangle(30);
    changePosition(gca,[0.06 0.03 0.02 -0.11]);
    put_axes_labels(gca,{[],[0 0 0]},{{'Correlation'},[0 0 0]});

    save_pdf(hf,mData.pdf_folder,'cell corr remap',600);
return;
end


%%
% out_CC = all_out_C.all_corr_an{3};
% out_AA = all_out_A.all_corr_an{1};
group = 1
if group == 1
    sel_out = mean_corr_C;
    xs = xs_C;
else
    sel_out = mean_corr_A;
    xs = xs_A;
end

ff = makeFigureRowsCols(107,[1 0.5 2 2],'RowsCols',[4 4],...
    'spaceRowsCols',[0.05 0.07],'rightUpShifts',[0.1 0.1],'widthHeightAdjustment',...
    [-95 -95]);
gg = 1;
set(gcf,'color','w');
set(gcf,'Position',[1 6 3.4 3.4]);
FS = mData.axes_font_size;
maskdisp = triu(ones(4,4),0);

for rr = 1:4
    for cc = 1:4
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
            set(gca,'YTick',[1 colsHalf cols],'YTickLabel',[ts(1)-2 ts(colsHalf) ts(cols)+2]);
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
save_pdf(ff.hf,mData.pdf_folder,sprintf('pos corr remap_%d',group),600);
return;
