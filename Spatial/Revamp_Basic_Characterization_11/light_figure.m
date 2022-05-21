function light_figure

while 1
    mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size;
    ei = evalin('base','ei'); 
    selContexts = [1 4 6]; rasterNames = {'light22T','light22T','light22T'};
    oL = get_data(ei,selContexts,rasterNames);
    Rs = oL.Rs;
%     resp_L = oL.resp;
    break;
end
n = 0;
si = [Lb ArL_L Lbs]
Rs = o.Rs(:,si);
%% light duration
for ii = 1:length(ei)
    tei = ei{ii};
    for jj = 1:length(tei.b.stim_r)
        dur_ligh(ii,jj) = tei.b.ts(tei.b.stim_f(jj))-tei.b.ts(tei.b.stim_r(jj));
    end
    for jj = 2:length(tei.b.stim_r)
        i_dur_ligh(ii,jj) = tei.b.ts(tei.b.stim_r(jj))-tei.b.ts(tei.b.stim_f(jj-1));
    end
end
%%
[resp_fractionC,resp_valsC,OIC,mean_OIC,resp_ORC,resp_OR_fractionC,resp_ANDC,resp_AND_fractionC] = get_responsive_fraction(Rs);
resp = get_cell_list(resp_valsC,[]);
out = find_population_vector_corr_remap(Rs,mR,resp_ORC);

n = 0;
%% Speed Figure

while 1
    for an = 1:size(Rs,1)
        mean_speed_over_trials(an,:) = nanmean(Rs{an,2}.speed);
    end
    hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1 1],'color','w');
    hold on; 
    xs = Rs{1,2}.xs; cis = Rs{1,2}.resp.cis;
    xs = xs - xs(cis(1,2));
    xticks = [1 19 38];
    mspeed = mean(mean_speed_over_trials(:,1:38)); semspeed = std(mean_speed_over_trials(:,1:38))/sqrt(5);
    plot(xs,mspeed);
    shadedErrorBar(xs,mspeed,semspeed);
    ylims = ylim;
    plot([xs(cis(1,2)) xs(cis(1,2))],[0 ylims(2)+3],'linewidth',0.1,'color','m');
    changePosition(gca,[0.13 0.15 -0.2 -0.35]);
    put_axes_labels(gca,{'Time (sec)',[0 0 0]},{{'Speed (cm/sec)'},[0 -4 0]});
    format_axes(gca);
    set(gca,'xtick',xs(xticks),'xticklabels',{'-2','0','2'});
    save_pdf(hf,mData.pdf_folder,'LED_speed',600);
    break;
end
%% dpca trying
while 1
    an = 1;
    xtck = [];
    t = [];
    for cn = 1:3
        raster = permute(Rs{an,cn}.sp_rasters1,[3 2 1]);
        xtck = [xtck raster(:,:)];
        if cn == 1
            t = [t Rs{an,cn}.xs];
        else
            t = [t Rs{an,cn}.xs+t(end)];
        end
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
    dpca__mine(Rs(an,:),psths);
    break;
end

%%
resp_valsC = resp_L.vals;
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
while 1
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
	break;
end
%%
Rs = o.Rs(:,[Lb_T Lbs_T]);
propsL = get_props_Rs(Rs,50);
an = 4; cn = 2;
% plotRasters_simplest(Rs{an,cn},find(propsL.good_FR{an,cn}))
% find(resp_valsC{an}(:,cn));
ff = makeFigureRowsCols(2020,[0.5 0.5 4 1],'RowsCols',[1 2],...
    'spaceRowsCols',[0.15 0.1],'rightUpShifts',[0.15 0.25],'widthHeightAdjustment',...
    [-200 -475]);
gg = 1;
set(gcf,'color','w');
set(gcf,'Position',[10 4 2 1]);
% ff = sample_rasters(Rs{an,cn},[558 328 168 55],ff);
ff = sample_rasters(Rs{an,cn},[80 269],ff);
set_obj(ff,{'xtick',[1 18 38],'xticklabels',[-2 0 2.2]})
save_pdf(ff.hf,mData.pdf_folder,sprintf('light_rastgers'),600);

%% for presentation
while 1
        %%
    si = [Lb_T Lbs_T];
    Rs = o.Rs(:,[Lb_T Lbs_T]);
    propsL = get_props_Rs(Rs,50);
    an = 4; cn = 1;
    % plotRasters_simplest(Rs{an,cn},find(propsL.good_FR{an,cn}))
    % find(resp_valsC{an}(:,cn));
    ff = makeFigureRowsCols(2020,[0.5 0.5 4 1],'RowsCols',[1 4],...
        'spaceRowsCols',[0.15 0.07],'rightUpShifts',[0.1 0.25],'widthHeightAdjustment',...
        [-100 -475]);
    gg = 1;
    set(gcf,'color','w');
    set(gcf,'Position',[10 4 4 1]);
    % ff = sample_rasters(Rs{an,cn},[558 328 168 55],ff);
%     ff = sample_rasters(Rs{an,cn},[148 59 100 110],ff);
    ff = sample_rasters(Rs{an,cn},[25 30 41 78],ff);
    set_obj(ff,{'xtick',[1 18 38],'xticklabels',[-2 0 2.2]})
    colormap parula
    save_pdf(ff.hf,mData.pdf_folder,sprintf('light_rasters_pre_inh'),600);
    
    
    %%
    break;
end

%% population vector and correlation single animal
an = 2;
si = si_light;
Rs = o.Rs(:,si);mR = o.mR(:,si);
ntrials = 50;
props1 = get_props_Rs(o.Rs,ntrials);
% si = si_light;
untuned = cell_list_op(props1.vals(:,si),[],'not');
tuned = props1.vals(:,si);
good_FR = cell_list_op(props1.good_FR(:,si),tuned,'and');
good_FR = props1.good_FR;
ff = makeFigureRowsCols(107,[1 0.5 4 1],'RowsCols',[2 3],...
    'spaceRowsCols',[0 -0.03],'rightUpShifts',[0.07 0.1],'widthHeightAdjustment',...
    [0.01 -60]);
set(gcf,'color','w');
set(gcf,'Position',[5 5 3.25 2]);
[CRc,aCRc,mRR] = find_population_vector_corr(Rs,mR,good_FR,0);
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
    out = out_C{1,1};
    ff = makeFigureRowsCols(106,[1 0.5 4 0.5],'RowsCols',[3 3],...
        'spaceRowsCols',[0.1 0.1],'rightUpShifts',[0.13 0.2],'widthHeightAdjustment',...
        [-150 -150]);
    set(gcf,'color','w');
    set(gcf,'Position',[5 3 5 5]);
    ff = show_remapping_corr_plots(mData,ff,out.mean_PV_corr,out.xs,[]);
    save_pdf(ff.hf,mData.pdf_folder,sprintf('remap_corr_C.pdf'),600);
end
%%
while 1
    good_FR = cell_list_op(props1.good_FR(:,si),props1.vals(:,si),'and');
    good_FR = props1.good_FR(:,si);
    resp_fractionC = exec_fun_on_cell_mat(good_FR,'sum')./exec_fun_on_cell_mat(good_FR,'length');
    dataT = array2table(resp_fractionC*100);
    dataT.Properties.VariableNames = {'L1','L2','L3'};
    within = array2table([1 2 3]');
    within.Properties.VariableNames = {'Cond'};
    within.Cond = categorical(within.Cond);
    ra = RMA(dataT,within);
    ra.ranova
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Cond','hsd'},[1.5 1 1]);
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
    break;
end

