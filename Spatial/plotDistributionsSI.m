function plotDistributionsSI(fn,data,nbins)
% close all
thisCols_all = {[0 0 0],'b','r',[0 0.7 0.3],'m','c'};

allSI = [];
for ii = 1:length(data)
    allSI = [allSI data{ii}.vals];
end
maxSI = max(allSI);
minSI = min(allSI);
bins = linspace(minSI,maxSI,nbins);
allBars = [];
for ii = 1:length(data)
    [bar1 xs] = hist(data{ii}.vals,bins); bar1 = 100*bar1/sum(bar1);
    allBars = [allBars;bar1];
end


ff = makeFigureWindow__one_axes_only(fn,[6 4 7 5],[0.13 0.27 0.85 0.7]);
axes(ff.ha);hold on;
hb = bar(xs,allBars');
for ii = 1:length(hb)
    set(hb(ii),'facecolor',thisCols_all{ii});
end
xlim([bins(1) bins(end)]);
set(gca,'TickDir','out','FontSize',14,'FontWeight','Bold');
xlabel('Values');
ylabel('Percentage');

% x1 = 5; x2 = x1+3; y1 = (60:-4.5:0); y1 = y1(1:4); y2 = y1;
% legendFontSize = 11;
% legs = {'nbb_air','nbb_belt','bc_air_belt','bb_air_belt'};
% for ii = 1:length(legs)
%     plot([x1 x2],[y1(ii) y2(ii)],'color',thisCols_all{ii},'linewidth',2);
%     text(x2+0.5,y1(ii),sprintf('%s',legs{ii}),'Color',thisCols_all{ii},'FontSize',legendFontSize);
% end
axesPos = ff.pos + [0.6 0.35 0 0];
axesPos(3:4) = [0.15 0.2];
axes('Position',axesPos);hold on;

for ii = 1:length(data)
    cs = cumsum(allBars(ii,:));
    plot(xs,cs,'color',thisCols_all{ii},'linewidth',1.5);
end
xlim([bins(1) bins(end)]);
ylim([0 100]);
set(gca,'TickDir','out','FontSize',8,'FontWeight','Bold');
% xlabel('MI(z-score)');
% ylabel('%');

    annovaVar = [];
    for ii = 1:length(data)
        annovaVar = [annovaVar data{ii}.vals'];
        names{ii} = data{ii}.name;
        [mVals(ii) semVals(ii)] = findMeanAndStandardError(data{ii}.vals);
    end
    if length(data) > 2
        gn = 1:length(data);
        [p,tbl,stats] = anova1(annovaVar,gn,'off');%,names);
        % [p,tbl,stats] = kruskalwallis(y2wa,subst,'on');
%         figure(2001);
        [c,~,~,gnames] = multcompare(stats,'CType','hsd','Display','off');
        pdf = c(:,6);
        hdf = pdf<0.05;
        selGroups = 1:length(data);
        combs = nchoosek(1:length(selGroups),2);

        for ii = 1:size(combs,1)
            inds = ismember(c(:,[1 2]),[selGroups(combs(ii,1)) selGroups(combs(ii,2))],'rows');
            prT(ii,1) = c(inds,6);
        end
        hrT = prT<0.05;
    else
        selGroups = 1:length(data);
        combs = nchoosek(1:length(selGroups),2);
        [hrT,prT,ci,stat] = ttest2(annovaVar(:,1),annovaVar(:,2))
        n = 0;
    end

    % bar graph from anova
    % ff = makeFigureWindow__one_axes_only(2,[1 4 1.75 2],[0.2 0.35 0.7 0.6]);
    ff = makeFigureWindow__one_axes_only(3,[2 4 5 5],[0.2 0.17 0.75 0.81]);
    axes(ff.ha);
    plotBarsWithSigLines(mVals,semVals,combs,[hrT prT],'colors',thisCols_all(selGroups),'ySpacingFactor',10);
    xlim([0.4 0.6+length(selGroups)]);
%     ylim([0 max(mVals)]);
    hyl = ylabel('Average Values');
    pos = get(hyl,'Position');pos = pos + [+0.1 0 0];set(hyl,'Position',pos);
    set(ff.ha,'linewidth',1.25);
    set(ff.ha,'TickDir','out','FontSize',9);
    set(ff.ha,'XTick',(1:length(data)));
    xtickangle(15);
