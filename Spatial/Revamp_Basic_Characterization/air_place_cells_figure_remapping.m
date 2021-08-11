function air_place_cells_figure_remapping

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

%% all cells
selRespA = get_cell_list(resp_valsC,[NaN NaN NaN]);
[all_corr_CA,all_corr_cell_CA,mean_corr_CA,mean_cell_corr_CA,xs_CA,paramA] = find_population_vector_corr_remap(Rs,mR,selRespA);
FDA = paramA.FD;
%% place cells in condition 1 only
selResp1 = get_cell_list(resp_valsC,[1 -2 -3]);
[all_corr_C1,all_corr_cell_C1,mean_corr_C1,mean_cell_corr_C1,xs_C1,param1] = find_population_vector_corr_remap(Rs,mR,selResp1);
FD1 = param1.FD;
%% place cells in condition 2 only
selResp2 = get_cell_list(resp_valsC,[-1 2 -3]);
[all_corr_C2,all_corr_cell_C2,mean_corr_C2,mean_cell_corr_C2,xs_C2,param2] = find_population_vector_corr_remap(Rs,mR,selResp2);
FD2 = param2.FD;
%% place cells in condition 3 only
selResp3 = get_cell_list(resp_valsC,[-1 -2 3]);
[all_corr_C3,all_corr_cell_C3,mean_corr_C3,mean_cell_corr_C3,xs_C3,param3] = find_population_vector_corr_remap(Rs,mR,selResp3);
FD3 = param3.FD;

n = 0;
%%
% varNames = {'Group','R11','R12','R13','R22','R23','R33'};
% 
% tab1 = ([squeeze(FD1(1,1,:)),squeeze(FD1(1,2,:)),squeeze(FD1(1,3,:)),squeeze(FD1(2,2,:)),squeeze(FD1(2,3,:)),squeeze(FD1(3,3,:))]);
% tab2 = ([squeeze(FD2(1,1,:)),squeeze(FD2(1,2,:)),squeeze(FD2(1,3,:)),squeeze(FD2(2,2,:)),squeeze(FD2(2,3,:)),squeeze(FD2(3,3,:))]);
if 1
% varNames = {'Group','R11','R22','R33'};
varNames = {'Group','C12','C13','C23'};

tabA = ([squeeze(FDA(1,1,:)),squeeze(FDA(2,2,:)),squeeze(FDA(3,3,:))]);
tab1 = ([squeeze(FD1(1,2,:)),squeeze(FD1(1,3,:)),squeeze(FD1(2,3,:))]);
tab2 = ([squeeze(FD2(1,2,:)),squeeze(FD2(1,3,:)),squeeze(FD2(2,3,:))]);
tab3 = ([squeeze(FD3(1,2,:)),squeeze(FD3(1,3,:)),squeeze(FD3(2,3,:))]);

tab1 = ([squeeze(mean_cell_corr_C1(1,2,:)),squeeze(mean_cell_corr_C1(1,3,:)),squeeze(mean_cell_corr_C1(2,3,:))]);
tab2 = ([squeeze(mean_cell_corr_C2(1,2,:)),squeeze(mean_cell_corr_C2(1,3,:)),squeeze(mean_cell_corr_C2(2,3,:))]);
tab3 = ([squeeze(mean_cell_corr_C3(1,2,:)),squeeze(mean_cell_corr_C3(1,3,:)),squeeze(mean_cell_corr_C3(2,3,:))]);



dataT = array2table([([ones(5,1);(2*ones(5,1));(3*ones(5,1))]) [tab1;tab2;tab3]]);
dataT.Properties.VariableNames = varNames;
dataT.Group = categorical(dataT.Group);
within = array2table([[1 2 3]']);
within.Properties.VariableNames = {'Cond'};
within.Cond = categorical(within.Cond);


ra = repeatedMeasuresAnova(dataT(6:10,2:end),within,0.05);

mVar = ra.est_marginal_means.Mean;semVar = ra.est_marginal_means.Formula_StdErr;
    combs = ra.mcs.combs; p = ra.mcs_wf.p; h = ra.mcs.p < 0.05;
%     row = [1 2]; ii = ismember(combs,row,'rows'); p(ii) = mcGroup{1,6}; h(ii) = 1; 
    % row = [5 6]; ii = ismember(combs,row,'rows'); p(ii) = mcTI{2,6}; h(ii) = 1; 
    % row = [3 4]; ii = ismember(combs,row,'rows'); p(ii) = mcTI{3,6}; h(ii) = 1; 

    xdata = [1:length(mVar)]; 
    colors = mData.colors;
    hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 3.44 1],'color','w');
    hold on;
%     tcolors = {colors{1};colors{1};colors{2};colors{2};colors{3};colors{3};colors{4};colors{4}};
    tcolors = repmat({colors{1};colors{2};colors{3}},4,1);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.25,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.001,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.7,'sigLinesStartYFactor',0.05);
%     for ii = 5:length(hbs)
%     set(hbs(ii),'facecolor','none','edgecolor',tcolors{ii});
%     end
    % plot([0.5 11],[-0.5 0.5],'linewidth',1.5)
%     maxY = maxY + 5 + 2;
    set(gca,'xlim',[0.25 xdata(end)+0.75],'ylim',[0 maxY],'FontSize',6,'FontWeight','Bold','TickDir','out');
    xticks = xdata(1:end); xticklabels = varNames(2:end);
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
%     xtickangle(30)
%     changePosition(gca,[0.2 0.03 -0.3 -0.05]);
    put_axes_labels(gca,{[],[0 0 0]},{{'Correlation'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('remap_FD'),600);
    return;
end
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

    sel_out = mean_corr_C1;
    xs = xs_C1;
    
    sel_out = param1.all_bw;
    xs = xs_C1;

ff = makeFigureRowsCols(101,[1 0.5 2 2],'RowsCols',[3 3],...
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
%         hc = putColorBar(ff.h_axes(rr,cc),[0.0 0.03 0 -0.05],[minC maxC],5,'eastoutside',[0.07 0.08 0.1 0.13]);
    end
end
save_pdf(ff.hf,mData.pdf_folder,sprintf('pos corr remap'),600);
return;
