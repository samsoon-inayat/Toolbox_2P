function figure1_Distributions
mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size;
pop_corr_C = load('pop_corr_C.mat');
pop_corr_A = load('pop_corr_A.mat');
tril_val = -20;
for ii = 1:4
    for an = 1:3
        corrmat = pop_corr_C.avg_C_conds{ii}(:,:,an);
        fcorrmat_C{an,ii} = imgaussfilt(corrmat,2);
%         wcp{an,ii} = get_width_corr_param(corrmat,fcorrmat_C{an,ii});
        maskmat = ones(size(corrmat));
        maskmattril = tril(maskmat,tril_val);% & ~tril(maskmat,-2);
        avgpc_C(an,ii) = nanmean(corrmat(maskmattril==1));
    end
end

for ii = 1:4
    for an = 1:5
        corrmat = pop_corr_A.avg_C_conds{ii}(:,:,an);
        fcorrmat_A{an,ii} = imgaussfilt(corrmat,2);
        maskmat = ones(size(corrmat));
        maskmattril = tril(maskmat,tril_val);% & ~tril(maskmat,-2);
        avgpc_A(an,ii) = nanmean(corrmat(maskmattril==1));
    end
end
n= 0;
%%
    data = avgpc_C;
    cmdTxt = sprintf('dataT_C = table(');
    for ii = 1:(size(data,2)-1)
        cmdTxt = sprintf('%sdata(:,%d),',cmdTxt,ii);
    end
    cmdTxt = sprintf('%sdata(:,size(data,2)));',cmdTxt);
    eval(cmdTxt);
    data = avgpc_A;
    cmdTxt = sprintf('dataT_A = table(');
    for ii = 1:(size(data,2)-1)
        cmdTxt = sprintf('%sdata(:,%d),',cmdTxt,ii);
    end
    cmdTxt = sprintf('%sdata(:,size(data,2)));',cmdTxt);
    eval(cmdTxt);
    dataT = [dataT_C;dataT_A]
    dataT.Properties.VariableNames = {'C1','C2','C3','C4'};
    dataT = [table([ones(3,1);2*ones(5,1)]) dataT];
    dataT.Properties.VariableNames{1} = 'Group';
    dataT.Group = categorical(dataT.Group)
    numCols = 4;
    colVar1 = [ones(1,numCols) 2*ones(1,numCols) 3*ones(1,numCols) 4*ones(1,numCols)];    colVar2 = [1:numCols 1:numCols 1:numCols 1:numCols];
    colVar1 = [1 2 3 4];
    within = table(colVar1');
%     within = table(colVar1',colVar2');
    within.Properties.VariableNames = {'Condition'};
%     within.Properties.VariableNames = {'Condition','Raster'};
    within.Condition = categorical(within.Condition);
%     within.Raster = categorical(within.Raster);
    ra = repeatedMeasuresAnova(dataT,within);

    mVar = ra.est_marginal_means.Mean;semVar = ra.est_marginal_means.Formula_StdErr;
    combs = ra.mcs.combs; p = ra.mcs.p; h = ra.mcs.p < 0.05;
     xdata = [1 2 3 4 6:9]; 
%     xdata = [1 2 3 4];
    colors = mData.colors;
    hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 2 1],'color','w');
    hold on;
    tcolors ={colors{1};colors{2};colors{3};colors{4};colors{1};colors{2};colors{3};colors{4}};
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.001,'sigTestName','','sigLineWidth',0.25,'BaseValue',0,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.7,'sigLinesStartYFactor',0);
    for ii = 5:length(hbs)
        set(hbs(ii),'facecolor','none','edgecolor',tcolors{ii});
    end
    % plot([0.5 11],[-0.5 0.5],'linewidth',1.5)
    set(gca,'xlim',[0.25 xdata(end)+0.75],'ylim',[-0.15 maxY],'FontSize',6,'FontWeight','Bold','TickDir','out');
    xticks = xdata(1:end)+0; xticklabels = {'C1-AD','C2-AD','C3-AD','C4-AD','C1-AD','C2-AD','C3-AD','C4-AD'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
    xtickangle(30);
    changePosition(gca,[0.075 0.0 0 -0.11]);
    put_axes_labels(gca,{[],[0 0 0]},{{'Percentage of','spatially tuned cells'},[0 -5 0]});
    
    save_pdf(hf,mData.pdf_folder,sprintf('Pop_vec_corr.pdf'),600);
    
function out = get_width_corr_param(pcm,fpcm)
n = 0;
threshold = 0.3;
for ii = 1:size(pcm,1)
    thisRow = fpcm(ii,:);
    [~,max_pos] = max(thisRow);
    thisRow_after_peak = thisRow(max_pos:end);
    pos = find(thisRow_after_peak<threshold,1,'first');
    if isempty(pos)
        wrp(ii) = NaN;
    else
        wrp(ii) = pos-max_pos;
    end
end
figure(1000);clf;
plot(fpcm)