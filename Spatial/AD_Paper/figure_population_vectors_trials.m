function figure_population_vectors_trials(fn,allRs,ccs)

mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size;
ei_C = evalin('base','ei10_C2'); 
ei_A = evalin('base','ei10_A2'); 

selContexts = [1 2 3 4];
rasterNames = {'airD','airD','airD','airD'};

RsC = get_rasters_data(ei_C,selContexts,rasterNames);
RsC = find_responsive_rasters(RsC,1:10);
[resp_fractionC,resp_valsC,OIC,mean_OIC,resp_ORC,resp_OR_fractionC,resp_ANDC,resp_AND_fractionC] = get_responsive_fraction(RsC);

CN = 1;
a_trials = {1:2,3:4,5:6,7:8,9:10};
mRsCt = [];
RsCt = [];
for ii = 1:length(a_trials)
    RsCt = [RsCt RsC(:,CN)];
    mRsCt = [mRsCt calc_mean_rasters(RsC(:,CN),a_trials{ii})];
end

[all_corr_C,all_corr_cell_C,mean_corr_C,mean_cell_corr_C,xs_C] = find_population_vector_corr_remap(RsCt,mRsCt,resp_ORC);

RsA = get_rasters_data(ei_A,selContexts,rasterNames);
RsA = find_responsive_rasters(RsA,1:10);
[resp_fractionA,resp_valsA,OIA,mean_OIA,resp_ORA,resp_OR_fractionA,resp_ANDA,resp_AND_fractionA] = get_responsive_fraction(RsA);

mRsAt = [];
RsAt = [];
for ii = 1:length(a_trials)
    RsAt = [RsAt RsA(:,CN)];
    mRsAt = [mRsAt calc_mean_rasters(RsA(:,CN),a_trials{ii})];
end

[all_corr_A,all_corr_cell_A,mean_corr_A,mean_cell_corr_A,xs_A] = find_population_vector_corr_remap(RsAt,mRsAt,resp_ORA);

n = 0;

%%
all_CC_C = mean_corr_C;
all_CC_A = mean_corr_A;
mean_cell_corr_C = mean_cell_corr_C;
mean_cell_corr_A = mean_cell_corr_A;
group = 2; 
%%
if group == 1
    sel_out = all_CC_C;
    sel_p.xs = xs_C;
else
    sel_out = all_CC_A;
    sel_p.xs = xs_A;
end

ff = makeFigureRowsCols(107,[1 0.5 2 2],'RowsCols',[5 5],...
    'spaceRowsCols',[0.01 0.01],'rightUpShifts',[0.15 0.12],'widthHeightAdjustment',...
    [-45 -45]);
gg = 1;
set(gcf,'color','w');
set(gcf,'Position',[1 6 1.7 1.7]);
FS = mData.axes_font_size;
maskdisp = triu(ones(5,5),0);

for rr = 1:5
    for cc = 1:5
        if ~maskdisp(rr,cc)
            delete(ff.h_axes(rr,cc));
            continue;
        end
        axes(ff.h_axes(rr,cc));
        this_pos_corr = sel_out{rr,cc};
        imagesc(this_pos_corr);
        box off;
        set(gca,'Ydir','Normal','linewidth',1,'FontSize',FS-1,'FontWeight','Bold');
        if rr == 1 && cc == 1
            cols = size(this_pos_corr,2);
            colsHalf = round(cols/2);
            ts = round(sel_p.xs);
            set(gca,'YTick',[1 colsHalf cols],'YTickLabel',[0 ts(colsHalf) ts(cols)+3]);
            h = ylabel('Pos (cm)');%    changePosition(h,[0 0 0]);
            set(gca,'XTick',[]);
        elseif rr == 5 && cc == 5
            h = xlabel('Pos (cm)');%    changePosition(h,[0 0 0]);
            cols = size(this_pos_corr,2);
            colsHalf = round(cols/2);
            ts = round(sel_p.xs);
            set(gca,'XTick',[1 colsHalf cols],'XTickLabel',[0 ts(colsHalf) ts(cols)+3]);
            set(gca,'YTick',[]);
        else
            axis off
        end
        minC = min(this_pos_corr(:));
        maxC = max(this_pos_corr(:));
%         text(5,cols-5,sprintf('(%.1f, %.1f)',minC,maxC),'FontSize',5,'Color','w');
%         hc = putColorBar(ff.h_axes(rr,cc),[0.0 0.03 0 -0.05],[minC maxC],5,'eastoutside',[0.07 0.11 0.1 0.16]);
    end
end
save_pdf(ff.hf,mData.pdf_folder,sprintf('pos corr trials_group_%d_cond_%d',group,CN),600);
%%
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
dataT.Properties.VariableNames = {'Group','T1234','T3456','T5678','T78910'};
dataT.Group = categorical(dataT.Group);
colVar1 = [1 2 3 4];    
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
xdata = [1 2 3 4 6:9];
colors = mData.colors;
hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.75 1],'color','w');
hold on;
tcolors ={colors{1};colors{2};colors{3};colors{4};colors{1};colors{2};colors{3};colors{4}};
[hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
    'ySpacing',0.1,'sigTestName','','sigLineWidth',0.25,'BaseValue',0,...
    'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.7,'sigLinesStartYFactor',0.1);
for ii = 5:length(hbs)
    set(hbs(ii),'facecolor','none','edgecolor',tcolors{ii});
end
% plot([0.5 11],[-0.5 0.5],'linewidth',1.5)
set(gca,'xlim',[0.25 xdata(end)+0.75],'ylim',[0 maxY],'FontSize',6,'FontWeight','Bold','TickDir','out');
xticks = xdata(1:end)+0; xticklabels = {'TC1','TC2','TC3','TC4','TC1','TC2','TC3','TC4'};
set(gca,'xtick',xticks,'xticklabels',xticklabels);
xtickangle(30);
changePosition(gca,[0.06 0.03 0.02 -0.11]);
put_axes_labels(gca,{[],[0 0 0]},{{'Correlation'},[0 0 0]});

save_pdf(hf,mData.pdf_folder,sprintf('cell corr remap trials _ group_%d_cond_%d',group,CN),600);
