function figure_place_remapping_AD(fn,allRs,ccs)

protocol_C = '10_C';
% protocol_A = '10_A';
ei_C = evalin('base','ei10_C');
% ei_A = evalin('base','ei10_A');
mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size;
ET_C = evalin('base',sprintf('ET%s',protocol_C));
% ET_A = evalin('base',sprintf('ET%s',protocol_A));
selAnimals_C = 1:length(ei_C)
% selAnimals_A = 1:length(ei_A)

% in the following variable all the measurements are in the matrices form
% for each variable colums indicate raster and stim marker types specified 
% the rows indicate condition numbers.
paramMs_C = parameter_matrices_ctrl('get','10_C_ctrl');
% paramMs_A = parameter_matrices('get','10_A');
% after getting all matrics, we can apply selection criteria to select a
% subgroup of cells
% here is the selection criteria in make_selC_structure function
cellsOrNot = 1; planeNumber = NaN; zMI_Th = 5; fwids = [1 120]; fcens = [0 140]; rs_th = 0.9;
% cellsOrNot = 1; planeNumber = NaN; zMI_Th = NaN; fwids = NaN; fcens = NaN; rs_th = NaN;
included_cells = 'All';
conditionsAndRasterTypes = [11;21;31;41];
selC = make_selC_struct(cellsOrNot,planeNumber,conditionsAndRasterTypes,zMI_Th,fwids,fcens,rs_th,NaN,NaN);
[cpMs_C,pMs_C] = parameter_matrices_ctrl('select','10_C',{paramMs_C,selC});
% [cpMs_A,pMs_A] = parameter_matrices('select','10_A',{paramMs_A,selC});
% perc_cells_C = parameter_matrices('print','10_C',{cpMs_C,pMs_C,ET_C,selAnimals_C});
% perc_cells_A = parameter_matrices('print','10_A',{cpMs_A,pMs_A,ET_A,selAnimals_A});
out_C = find_corr_coeff_rasters(pMs_C',paramMs_C,selAnimals_C,ei_C,conditionsAndRasterTypes',selC,cpMs_C);
% out_A = find_corr_coeff_rasters(pMs_A',paramMs_A,selAnimals_A,ei_A,conditionsAndRasterTypes',selC,cpMs_A);
n = 0;

%%
if 1
    data_C_motion = out_C.mean_over_cells;
    data_A_motion = out_A.mean_over_cells;
    dataT_C = array2table(data_C_motion); dataT_C.Properties.VariableNames = {'CC1','CC2','CC3'};
    dataT_A = array2table(data_A_motion); dataT_A.Properties.VariableNames = {'CC1','CC2','CC3'};
    groupT = table([ones(size(dataT_C,1),1);2*ones(size(dataT_A,1),1)]); groupT.Properties.VariableNames = {'Group'};
    between = [groupT [dataT_C;dataT_A]]; between.Group = categorical(between.Group);
    within = table([1 2 3]');
    within.Properties.VariableNames = {'dCond'};
    within.dCond = categorical(within.dCond);
    rmaR = repeatedMeasuresAnova(between,within);
    writetable(between,fullfile(mData.pdf_folder,sprintf('%s_dataT_%s.xlsx','figure_place_remapping_AD',included_cells)));
    writetable(rmaR.ranova,fullfile(mData.pdf_folder,sprintf('%s_dataT_stats_output_%s.xlsx','figure_place_remapping_AD',included_cells)),'WriteRowNames',true);
    
%     findMeanAndStandardError
    mVar = rmaR.est_marginal_means.Mean;
    semVar = rmaR.est_marginal_means.Formula_StdErr;
    combs = rmaR.combs;
    p = rmaR.p; h = p<0.05;
    xdata = [1 2 3 4 5 6]; maxY = 0.3;
    colors = mData.colors;
    hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 2.25 1],'color','w');
    hold on;
    tcolors = {colors{1};colors{2};colors{3};colors{1};colors{2};colors{3}};
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'maxY',maxY,'ySpacing',0.1,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.7,'sigLinesStartYFactor',0.1);
    for ii = 4:length(hbs)
        set(hbs(ii),'facecolor','none','edgecolor',tcolors{ii});
    end
    % plot([0.5 11],[-0.5 0.5],'linewidth',1.5)
    set(gca,'xlim',[0.25 xdata(end)+0.75],'ylim',[0 maxY+0],'FontSize',6,'FontWeight','Bold','TickDir','out');
    xticks = xdata(1:end); xticklabels = {'Ctrl-12','Ctrl-23','Ctrl-34','AD-12','AD-23','AD-34'};xtickangle(45)
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
    changePosition(gca,[0.02 0.03 0.02 -0.11]);
    put_axes_labels(gca,{[],[0 0 0]},{'Correlation Coeff.',[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('%s_bargraph_overall_%s','figure_place_remapping_AD',included_cells),600);
end
n = 0;

%%
runthis = 1;
if runthis
    cond = 3;
    for ii = 1:length(out_C.all_css_C)
        data_C(ii) = {squeeze(out_C.all_css_C{ii}(cond,:))};
    end
    for ii = 1:length(out_A.all_css_C)
        data_A(ii) = {squeeze(out_A.all_css_C{ii}(cond,:))};
    end
    data = [data_C' data_A'];
    minBin = min([out_C.gAllVals out_A.gAllVals]);
    maxBin = max([out_C.gAllVals out_A.gAllVals]);
    incr = (maxBin-minBin)/100;
    hf = figure(1000);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 6 2 1.5],'color','w');
    hold on;
    [ha,hb,hca,sigR] = plotDistributions(data,'colors',colors,'maxY',90,'cumPos',[0.5 0.26 0.25 0.5],'min',minBin,'incr',incr,'max',maxBin);
    hold on;
    legs = {sprintf('%s-C','CC'),sprintf('%s-A','CC')};
%     ylim([0 10]);
%     xlim([-3 12])
    xlims = xlim; dx = xlims(2) - xlims(1); ylims = ylim; dy = ylims(2) - ylims(1);
    legs{length(legs)+1} = [xlims(1)+dx/1.5 dx/30 ylims(1)+dy/3 dy/15];
    putLegend(ha,legs,'colors',colors);
    axes(ha);
    h = xlabel('Pearson Correlation Coeff.');%changePosition(h,[0 -dy/3 0]);
    h = ylabel('Percentage');changePosition(h,[-0.1 0 0]);
    set(gca,'FontSize',axes_font_size,'FontWeight','Bold');changePosition(ha,[0.07 0.01 -0.05 -0.05]);
    file_name = sprintf('%s_distribution CC_%s_%d','figure_place_remapping_AD',included_cells,cond);
    save_pdf(hf,mData.pdf_folder,file_name,600);
return;
end



