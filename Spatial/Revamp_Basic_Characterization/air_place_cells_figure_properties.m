function air_place_cells_figure_properties

mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size;
ei_11_15 = evalin('base','ei_11_15'); 
ei_2_3 = evalin('base','ei_2_3'); 

selContexts = [3 4 5];
rasterNames = {'airD','airD','airD'};
% Rs = get_rasters_data(ei_11_15,selContexts,rasterNames);
Rs = get_rasters_data(ei_2_3,selContexts,rasterNames);
Rs = find_responsive_rasters(Rs,1:10);
mR = calc_mean_rasters(Rs,1:10);
[CRc,aCRc,mRR] = find_population_vector_corr(Rs,mR,1);
[resp_fractionC,resp_valsC,OIC,mean_OIC,resp_ORC,resp_OR_fractionC,resp_ANDC,resp_AND_fractionC] = get_responsive_fraction(Rs);
[all_corr_C,all_corr_cell_C,mean_corr_C,mean_cell_corr_C,xs_C] = find_population_vector_corr_remap(Rs,mR,resp_ORC);

n = 0;

%%

all_variables = {'all_zMIs','all_fFR','all_fwidths','all_frs',''};
ylabels = {{'Mutual Information','(z-score)'},{'Peak Firing','Rate (AU)'},'Field Widths (cm)','R-Squared',{'Spatially Tuned', 'Cells (%)'}};

vn = 3;

number_of_bins = 4;
[all_valsC,all_vals_NC] = get_values(Rs,number_of_bins,all_variables{vn});

if 0
    all_valsC = all_vals_NC;
    vn = 5;
end
varNames = [];
if number_of_bins > 1
    ind = 1;
    w1 = [];
    w2 = [];
    for ii = 1:(size(all_valsC,2)/number_of_bins)
        for jj = 1:number_of_bins
            varNames{ind} = sprintf('C%dB%d',ii,jj);
            xticklabels{ind} = sprintf('B%d',jj);
            temp_tcolors{ind} = colors{ii};
            w1 = [w1 ii];
            w2 = [w2 jj];
            ind = ind + 1;
        end
    end

    dataT = array2table([all_valsC]);
    dataT.Properties.VariableNames = {varNames{:}};
    within = array2table([w1' w2']);
    within.Properties.VariableNames = {'Cond','Bin'};
    within.Cond = categorical(within.Cond);
    within.Bin = categorical(within.Bin);
    ra = repeatedMeasuresAnova(dataT,within,0.05);
    ra.ranova
else
    temp_tcolors = repmat(colors(1:4),2,1);
    varNames = {'C3','C4','C5'};
    xticklabels = varNames;
    dataT = array2table(all_valsC);
    dataT.Properties.VariableNames = {varNames{:}};
    within = array2table([1 2 3]');
    within.Properties.VariableNames = {'Cond'};
    within.Cond = categorical(within.Cond);
    ra = repeatedMeasuresAnova(dataT,within,0.05);
    ra.ranova
end
% writetable(dataT,fullfile(mData.pdf_folder,sprintf('%s_values.xlsx',all_variables{vn})));

n = 0;
% average distributions w.r.t centers for the two groups
if 1
    mVar = ra.est_marginal_means.Mean; semVar = ra.est_marginal_means.Formula_StdErr;
    tcolors = [temp_tcolors temp_tcolors];
    combs = ra.mcs.combs; p = ra.mcs.p; h = p<0.05;
    xdata = [1:length(mVar)];
    colors = mData.colors;
%     hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 6.9 2],'color','w');
    hf = figure(6);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 5 3.5 1],'color','w');
    hold on;
    [hbs,maxY]  = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.1,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.5,'sigLinesStartYFactor',0.01);

    % plot([0.5 11],[-0.5 0.5],'linewidth',1.5)
    set(gca,'xlim',[0.25 xdata(end)+0.15],'ylim',[0 maxY+1],'FontSize',6,'FontWeight','Bold','TickDir','out');
    xticks = xdata;%(1:2:end)+0.5; 
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
    changePosition(gca,[0.05 0.03 -0.1 -0.11]);
    put_axes_labels(gca,{[],[0 0 0]},{ylabels{vn},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('%s_distributions_over_belt_%d',all_variables{vn},number_of_bins),600);
    return;
end
%%
if 1
%%
mVar = ra.est_marginal_means_wf.Mean;
semVar = ra.est_marginal_means_wf.Formula_StdErr;
combs = ra.mcs_wf.combs; p = ra.mcs_wf.p; h = ra.mcs_wf.p < 0.05;
xdata = [1:length(mVar)];
tcolors = [temp_tcolors];
hf = figure(6);clf;set(gcf,'Units','Inches');set(gcf,'Position',[10 7 3.25 1],'color','w');
hold on;
% tcolors = repmat(tcolors,2,1)';
[hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
    'ySpacing',0.01,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.001,...
    'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.01);

set(gca,'xlim',[0.25 xdata(end)+0.75],'ylim',[0 maxY],'FontSize',6,'FontWeight','Bold','TickDir','out','xcolor','k','ycolor','k');
xticks = xdata; xticklabels = {'B1','B2','B3'}; xticklabels = repmat(xticklabels,1,2);
set(gca,'xtick',xticks,'xticklabels',xticklabels);

changePosition(gca,[0.0 0.02 0 -0.011])
put_axes_labels(gca,{[],[0 0 0]},{'',[0 0 0]});

save_pdf(hf,mData.pdf_folder,sprintf('%s_distributions_over_belt_%d_pooled',all_variables{vn},number_of_bins),600);
end

function [all_vals,all_vals_N] = get_values(Rs,number_of_bins,var)

all_vals = [];
all_vals_N = [];

for cc = 1:size(Rs,2)
    these_vals = NaN(size(Rs,1),number_of_bins);
    these_vals_N = NaN(size(Rs,1),number_of_bins);
    for rr = 1:size(Rs,1)
        if cc == 4
            n = 0;
        end
        R = Rs{rr,cc};
        mbl = mean(R.beltLength);
        binSize = mbl/number_of_bins;
        resp = R.resp.vals';
        [rs,as,bs,cs] = get_gauss_fit_parameters(R.gauss_fit_on_mean,R.bin_width);
        if strcmp(var,'all_zMIs')
            vals = R.info_metrics.ShannonMI_Zsh';
        end
        if strcmp(var,'all_fFR')
            vals = as';
            resp(vals>5000) = 0;
        end
        if strcmp(var,'all_fwidths')
            vals = cs';
        end
        if strcmp(var,'all_frs')
            vals = rs';
        end
        vals = vals(resp);
        if number_of_bins > 1
        bs = bs(resp);
        bins = 0:binSize:mbl;
        [N,E,Bi] = histcounts(bs,bins);
        mean_sv_Vals = [];
        for bb = 1:length(N)
            mean_sv_Vals(bb) = nanmean(vals(Bi == bb));
        end
        these_vals(rr,:) = mean_sv_Vals;
        these_vals_N(rr,:) = 100*(N/sum(N));
        else
            these_vals(rr,:) = nanmean(vals);
            these_vals_N(rr,:) = 100*length(vals)/length(resp);
        end
    end
    all_vals = [all_vals these_vals];
    all_vals_N = [all_vals_N these_vals_N];
end
