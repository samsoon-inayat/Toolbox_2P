function light_figure

mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size;
ei_11_15 = evalin('base','ei'); 

selContexts = [1 4 6];
rasterNames = {'light22T','light22T','light22T'};
Rs = get_rasters_data(ei_11_15,selContexts,rasterNames);
% Rs = get_rasters_data(ei_2_3,selContexts,rasterNames);
Rs = find_responsive_rasters(Rs,1:10);
mR = calc_mean_rasters(Rs,1:10);
[CRc,aCRc,mRR] = find_population_vector_corr(Rs,mR,1,0);
n = 0;
%%
[resp_fractionC,resp_valsC,OIC,mean_OIC,resp_ORC,resp_OR_fractionC,resp_ANDC,resp_AND_fractionC] = get_responsive_fraction(Rs);
resp = get_cell_list(resp_valsC,[]);
out = find_population_vector_corr_remap(Rs,mR,resp_ORC);

n = 0;
%% Speed Figure
if 1
    for an = 1:size(Rs,1)
        mean_speed_over_trials(an,:) = nanmean(Rs{an,2}.speed);
    end
    hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.5 1],'color','w');
    hold on; 
    xs = Rs{1,2}.xs; cis = Rs{1,2}.resp.cis;
    xs = xs - xs(cis(1,2));
    xticks = [cis(1,:) 38];
    mspeed = mean(mean_speed_over_trials(:,1:38)); semspeed = std(mean_speed_over_trials(:,1:38))/sqrt(5);
    plot(xs,mspeed);
    shadedErrorBar(xs,mspeed,semspeed);
    ylims = ylim;
    plot([xs(cis(1,2)) xs(cis(1,2))],[0 ylims(2)+3],'linewidth',0.1,'color','m');
    changePosition(gca,[0.07 0.15 -0.1 -0.15]);
    put_axes_labels(gca,{'Time (sec)',[0 0 0]},{{'Speed (cm/sec)'},[0 -4 0]});
    save_pdf(hf,mData.pdf_folder,'LED_speed',600);
end
%%
if 1
    an = 1;
    xtck = [];
    for cn = 1:3
        raster = permute(Rs{an,cn}.sp_rasters1,[3 2 1]);
        xtck = [xtck raster(:,:)];
    end
    N = size(xtck,1);
    T = size(Rs{an,1}.sp_rasters1,2);
    E = size(Rs{an,1}.sp_rasters1,1);
    C = 3;
    firing_rates = reshape(xtck,[N T E C]); firing_rates = permute(firing_rates,[1 4 2 3]);
    x = repmat(nanmean(xtck,2),[1 size(xtck,2)]);
    xt = reshape(xtck - x,[N C T E]); xt = nanmean(xt,3); xt = nanmean(xt,4); xt = repmat(xt,[1 1 T E]); xt = xt(:,:);
    xc = reshape(xtck - x,[N C T E]); xc = nanmean(xc,2); xc = nanmean(xc,4); xc = repmat(xc,[1 C 1 E]); xc = xc(:,:);
    xtc = reshape(xtck - x - xt - xc,[N C T E]); xtc = nanmean(xtc,4); xtc = repmat(xtc,[1 1 1 E]); xtc = xtc(:,:);
    psths = nanmean(firing_rates,4); psths = repmat(psths,[1 1 1 E]);
    psths = psths(:,:);
    etck = xtck - psths;
end

dpca__mine(Rs(an,:),psths);


%%

trials = mat2cell([1:10]',ones(size([1:10]')));
resp = get_cell_list(resp_valsC,[1;2;3]);
out1 = find_population_vector_corr_remap_trials(Rs(:,1),resp,trials);
out2 = find_population_vector_corr_remap_trials(Rs(:,2),resp,trials);
out3 = find_population_vector_corr_remap_trials(Rs(:,3),resp,trials);


n = 0;

%%

var1 = out1.adj_SP_corr_diag;
var2 = out2.adj_SP_corr_diag;
var3 = out3.adj_SP_corr_diag;
%%
if 1
    for rr = 1:size(var1,1)
        for cc = 1:size(var1,2)
            var_1(rr,cc) = nanmean(var1{rr,cc});
            var_2(rr,cc) = nanmean(var2{rr,cc});
            var_3(rr,cc) = nanmean(var3{rr,cc});
        end
    end
    ind = 1; ind_val = 1;
    for ii = 1:3
        for cc = 1:size(var_1,2)
            varNames{cc+(9*(ii-1))} = sprintf('C%d_T%d%d',ii,cc,cc+1);
            xlabels{cc+(9*(ii-1))} = sprintf('C%d-T%d-T%d',ii,cc,cc+1);
            xdata(ind) = ind_val;
            ind = ind+1; ind_val = ind_val + 1;
        end
        ind_val = ind_val+3;
    end
    
    dataT = array2table([[var_1 var_2 var_3]]);
    dataT.Properties.VariableNames = varNames;
    colVar1 = [1:size(var_1,2)]; colVar2 = [ones(size(colVar1)) 2*ones(size(colVar1)) 3*ones(size(colVar1))];
    colVar1 = repmat(colVar1,1,3);
    
    within = array2table([colVar2' colVar1']);
    within.Properties.VariableNames = {'Condition','TrialPairs'};
    within.TrialPairs = categorical(within.TrialPairs);
    within.Condition = categorical(within.Condition);
    ra = repeatedMeasuresAnova(dataT,within);
%     writetable(dataT,fullfile(mData.pdf_folder,'Remapping_Trials.xls'));

    mVar = ra.est_marginal_means.Mean;semVar = ra.est_marginal_means.Formula_StdErr;
    combs = ra.mcs.combs; p = ra.mcs.p; h = ra.mcs.p < 0.05;
    nbars = length(mVar)/3;
    colors = mData.colors;
    hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 5 1.25],'color','w');
    hold on;
    tcolors = colors(1:nbars); tcolors = repmat(tcolors,1,3);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.05,'sigTestName','','sigLineWidth',0.25,'BaseValue',0,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.7,'sigLinesStartYFactor',0.1);
%     for ii = (nbars+1):length(hbs)
%         set(hbs(ii),'facecolor','none','edgecolor',tcolors{ii});
%     end
    % plot([0.5 11],[-0.5 0.5],'linewidth',1.5)
    set(gca,'xlim',[0.25 xdata(end)+0.75],'ylim',[0 maxY],'FontSize',6,'FontWeight','Bold','TickDir','out');
    xticks = xdata(1:end)+0; xticklabels = xlabels; xticklabels = repmat(xticklabels,1,2);
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
    xtickangle(30);
    changePosition(gca,[-0.06 0.03 0.15 -0.05]);
    put_axes_labels(gca,{[],[0 0 0]},{{'Correlation'},[0 0 0]});

    save_pdf(hf,mData.pdf_folder,'remap bar graph_trials_light',600);
return;
end
%%
an = 1; cn = 2;
% plotRasters_simplest(Rs{an,cn})
% find(resp_valsC{an}(:,cn));
ff = makeFigureRowsCols(2020,[0.5 0.5 4 1],'RowsCols',[1 4],...
    'spaceRowsCols',[0.15 0.03],'rightUpShifts',[0.08 0.25],'widthHeightAdjustment',...
    [-50 -375]);
gg = 1;
set(gcf,'color','w');
set(gcf,'Position',[10 4 3.25 1]);
ff = sample_rasters(Rs{an,cn},[558 328 168 55],ff);
save_pdf(ff.hf,mData.pdf_folder,sprintf('light_rastgers'),600);
%% population vector and correlation single animal
an = 2;
ff = makeFigureRowsCols(107,[1 0.5 4 1],'RowsCols',[2 3],...
    'spaceRowsCols',[0 -0.03],'rightUpShifts',[0.07 0.1],'widthHeightAdjustment',...
    [0.01 -60]);
set(gcf,'color','w');
set(gcf,'Position',[5 5 3.25 2]);
ff = show_population_vector_and_corr(mData,ff,Rs(an,:),mRR(an,:),CRc(an,:),[],[]);
save_pdf(ff.hf,mData.pdf_folder,sprintf('population_vector_corr.pdf'),600);

%% average correlation of all animals
ff = makeFigureRowsCols(107,[1 0.5 4 0.5],'RowsCols',[1 3],...
    'spaceRowsCols',[0 -0.03],'rightUpShifts',[0.075 0.2],'widthHeightAdjustment',...
    [0.01 -220]);
set(gcf,'color','w');
set(gcf,'Position',[5 5 3.25 1]);
ff = show_population_vector_and_corr(mData,ff,Rs(an,:),[],aCRc,[],[]);
save_pdf(ff.hf,mData.pdf_folder,sprintf('average_population_vector_corr.pdf'),600);

%%
%% average correlation of all animals
if 1
    ff = makeFigureRowsCols(106,[1 0.5 4 0.5],'RowsCols',[3 3],...
        'spaceRowsCols',[0.1 0.1],'rightUpShifts',[0.13 0.2],'widthHeightAdjustment',...
        [-150 -150]);
    set(gcf,'color','w');
    set(gcf,'Position',[5 3 5 5]);
    ff = show_remapping_corr_plots(mData,ff,out.mean_PV_corr,out.xs,[]);
    save_pdf(ff.hf,mData.pdf_folder,sprintf('remap_corr_C.pdf'),600);
end
%%
dataT = array2table(resp_fractionC*100);
dataT.Properties.VariableNames = {'L1','L2','L3'};
within = array2table([1 2 3]');
within.Properties.VariableNames = {'Cond'};
within.Cond = categorical(within.Cond);
ra = repeatedMeasuresAnova(dataT,within,0.05);

%%
mVar = ra.est_marginal_means.Mean;semVar = ra.est_marginal_means.Formula_StdErr;
    combs = ra.mcs.combs; p = ra.mcs.p; h = ra.mcs.p < 0.05;
%     row = [1 2]; ii = ismember(combs,row,'rows'); p(ii) = mcGroup{1,6}; h(ii) = 1; 
    % row = [5 6]; ii = ismember(combs,row,'rows'); p(ii) = mcTI{2,6}; h(ii) = 1; 
    % row = [3 4]; ii = ismember(combs,row,'rows'); p(ii) = mcTI{3,6}; h(ii) = 1; 

    xdata = [1 2 3]; 
    colors = mData.colors;
    hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.25 1],'color','w');
    hold on;
%     tcolors = {colors{1};colors{1};colors{2};colors{2};colors{3};colors{3};colors{4};colors{4}};
    tcolors = {colors{1};colors{2};colors{3};colors{4};colors{1};colors{2};colors{3};colors{4};};
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',2,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.001,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.7,'sigLinesStartYFactor',0.05);
    for ii = 5:length(hbs)
    set(hbs(ii),'facecolor','none','edgecolor',tcolors{ii});
    end
    % plot([0.5 11],[-0.5 0.5],'linewidth',1.5)
%     maxY = maxY + 5 + 2;
    set(gca,'xlim',[0.25 xdata(end)+0.75],'ylim',[0 maxY],'FontSize',6,'FontWeight','Bold','TickDir','out');
    xticks = xdata(1:end); xticklabels = {'C1','C4','C1'''};
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
%     xtickangle(30)
    changePosition(gca,[0.2 0.03 -0.3 -0.05]);
    put_axes_labels(gca,{[],[0 0 0]},{{'Light Responsive','Cells (%)'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('fraction_light_responsive'),600);
%% overlap of cells
hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[7 7 1.25 1],'color','w');
    hold on;
    imagesc(mean_OIC);
    axis equal
    colorbar;
    xlim([0.5 3.5]);
    ylim([0.5 3.5]);
    changePosition(gca,[0.1 0.03 0 0]);
    xticklabels = {'C1','C4','C1'''};
    set(gca,'XTick',[1 2 3],'XTickLabels',xticklabels,'YTick',[1 2 3],'YTickLabels',(xticklabels));
    set(gca,'Ydir','reverse','linewidth',0.5,'FontSize',6,'FontWeight','Bold');
    save_pdf(hf,mData.pdf_folder,sprintf('light_overlap_image'),600);
%% overlap RM bar graph
hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.25 1],'color','w');
    hold on;
    for ii = 1:5
        C12(ii) = OIC{ii}(1,2);
        C13(ii) = OIC{ii}(1,3);
        C23(ii) = OIC{ii}(2,3);
    end
    dataT = array2table([C12' C13' C23']);
    dataT.Properties.VariableNames = {'L1','L2','L3'};
    within = array2table([1 2 3]');
    within.Properties.VariableNames = {'Cond'};
    within.Cond = categorical(within.Cond);
    ra = repeatedMeasuresAnova(dataT,within,0.05);

%%
mVar = ra.est_marginal_means.Mean;semVar = ra.est_marginal_means.Formula_StdErr;
    combs = ra.mcs.combs; p = ra.mcs.p; h = ra.mcs.p < 0.05;
%     row = [1 2]; ii = ismember(combs,row,'rows'); p(ii) = mcGroup{1,6}; h(ii) = 1; 
    % row = [5 6]; ii = ismember(combs,row,'rows'); p(ii) = mcTI{2,6}; h(ii) = 1; 
    % row = [3 4]; ii = ismember(combs,row,'rows'); p(ii) = mcTI{3,6}; h(ii) = 1; 

    xdata = [1 2 3]; 
    colors = mData.colors;
    hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.25 1],'color','w');
    hold on;
%     tcolors = {colors{1};colors{1};colors{2};colors{2};colors{3};colors{3};colors{4};colors{4}};
    tcolors = {colors{1};colors{2};colors{3};colors{4};colors{1};colors{2};colors{3};colors{4};};
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',2,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.001,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.7,'sigLinesStartYFactor',0.05);
    for ii = 5:length(hbs)
    set(hbs(ii),'facecolor','none','edgecolor',tcolors{ii});
    end
    % plot([0.5 11],[-0.5 0.5],'linewidth',1.5)
%     maxY = maxY + 5 + 2;
    set(gca,'xlim',[0.25 xdata(end)+0.75],'ylim',[0 maxY],'FontSize',6,'FontWeight','Bold','TickDir','out');
    xticks = xdata(1:end); xticklabels = {'C1-C4','C1-C1''','C4-C1'''};
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
    xtickangle(30)
    changePosition(gca,[0.2 0.03 -0.3 -0.05]);
    put_axes_labels(gca,{[],[0 0 0]},{{'Overlap Index'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('overlap_stats'),600);
    
