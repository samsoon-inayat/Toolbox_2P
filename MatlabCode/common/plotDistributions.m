function [ha,hb,hca,varargout] = plotDistributions (distD,varargin)
hb = NaN; hca = NaN; ha = NaN;
if nargin == 1
    if ismatrix(distD)
        distDo = [];
        allVals = [];
        allValsG = [];
        for cc = 1:size(distD,2)
            allValsGt = [];
            for rr = 1:size(distD,1)
                thisVal = distD{rr,cc};
                thisVal = nanmean(thisVal,2);
%                 thisVal = max(thisVal,[],2);
                thisVal = thisVal(:);
%                 thisVal = thisVal(thisVal > 0);
                allVals = [allVals;thisVal(:)];
                distDo{rr,cc} = thisVal;
                allValsGt = [allValsGt;thisVal(:)];
            end
            allValsG{cc} = allValsGt;
        end
        ha = distDo;
        hb = allVals;
        hca = allValsG;
    end
    return;
end



% allD = [];
% for ii = 1:length(distD)
%     allD = [allD distD{ii}];
% end
% maxB = ceil(max(allD));
% minB = floor(min(allD));

p = inputParser;
default_colors = distinguishable_colors(20);
default_ySpacingFactor = 10;
addRequired(p,'distD',@iscell);
addOptional(p,'incr',inf,@isnumeric);
addOptional(p,'min',-inf,@isnumeric);
addOptional(p,'max',inf,@isnumeric);
addOptional(p,'colors',default_colors,@iscell);
addOptional(p,'maxY',100,@isnumeric);
addOptional(p,'cumPos',0,@isnumeric);
addOptional(p,'barGraph',{},@iscell);
addOptional(p,'legend',{},@iscell);
addOptional(p,'BaseValue',0.2,@isnumeric);
addOptional(p,'do_mean','Yes');
parse(p,distD,varargin{:});

cols = p.Results.colors;
maxY = p.Results.maxY;
cumPos = p.Results.cumPos;
barGraph = p.Results.barGraph;
legs = p.Results.legend;
bv = p.Results.BaseValue;
do_mean = p.Results.do_mean;
% legs = temp(1:(end-1));
% specs = temp{end};

if p.Results.min ~= -inf
    minB = p.Results.min;
end

if p.Results.max ~= inf
    maxB = p.Results.max;
end

if p.Results.incr ~= inf
    incr = p.Results.incr;
else
    incr = ((maxB - minB)/20);
end

bins = (minB+incr):incr:(maxB-incr);
% else
%     bins = incr:incr:(maxB-incr);
% end
if isvector(distD) & strcmp(do_mean,'No')
    allBars = [];
    for ii = 1:length(distD)
        bd = distD{ii};
        [bar1 xs] = hist(bd,bins); bar1 = 100*bar1/sum(bar1);
        allBars = [allBars;bar1];
    end
    hb = bar(xs,allBars');
    for ii = 1:length(hb)
        set(hb(ii),'facecolor',cols{ii},'barwidth',0.7,'EdgeColor','none','BaseValue',bv,'ShowBaseline','off');
    end
    xlim([bins(1)-(incr/2) bins(end)+(incr/2)]);
    ylim([0 maxY]);
    set(gca,'TickDir','out','FontSize',7,'FontWeight','Normal','XTick',bins(1:2:end),'linewidth',0.25);
    ha = gca;
    % putLegend(gca,legs,specs,'colors',cols);
    if length(cumPos) > 1
        pos = get(gca,'Position');
        axesPos = pos + [cumPos(1) cumPos(2) 0 0];
        axesPos(3:4) = [cumPos(3) cumPos(4)];
        hca = axes('Position',axesPos);hold on;
        xlims = xlim;
        for ii = 1:length(distD)
            cs = cumsum(allBars(ii,:));
            plot(xs,cs,'color',cols{ii},'linewidth',0.5);
        end
        xlim([bins(1)-(incr/2) bins(end)+(incr/2)]);
        ylim([0 100]);
        set(gca,'TickDir','out','FontSize',5,'FontWeight','Normal');
        xlims = xlim; dx = diff(xlims);
        text(xlims(1)+dx/10,120,'Cumulative','FontSize',9,'FontWeight','Normal');
    end

    sigR = significanceTesting(distD);

    if ~isempty(barGraph)
        barPos = barGraph{2};
        maxY = barGraph{4};
        ySpacing = barGraph{6};
        sigTestName = barGraph{8};
        sigColor = barGraph{10};
        baseValue = barGraph{12};
        pos = get(gca,'Position');
        axesPos = pos + [barPos(1) barPos(2) 0 0];
        axesPos(3:4) = [barPos(3) barPos(4)];
        hca = axes('Position',axesPos);hold on;
        hs = sigR.anova.multcompare.h; ps = sigR.anova.multcompare.p;
        plotBarsWithSigLines(sigR.means,sigR.sems,sigR.combs,[hs ps],'colors',cols,'sigColor',sigColor,'maxY',maxY,'ySpacing',ySpacing,'sigTestName',sigTestName,'sigLineWidth',0.25,'BaseValue',baseValue);
        xlim([0.4 0.6+length(sigR.means)]);
        set(gca,'TickDir','out','FontSize',5,'FontWeight','Normal');
    end

    if nargout == 4
        varargout{1} = sigR;
    end

    axes(ha);
    if ~isempty(legs)
        putLegend(gca,legs,'colors',cols);
    end
end

if isvector(distD) & strcmp(do_mean,'Yes')
    allBars = [];
    for ii = 1:length(distD)
        bd = distD{ii};
        [bar1 xs] = hist(bd,bins); bar1 = 100*bar1/sum(bar1);
        allBars = [allBars;bar1];
    end
    [mVals,semVals] = findMeanAndStandardError(allBars);
    if length(mVals) > 1
%         shadedErrorBar(xs,mVals,semVals);
        hb = bar(xs,mVals);
        errorbar(xs,mVals,semVals,'.');
    else
        hb = bar(xs,allBars');
    end
%     for ii = 1:length(hb)
%         set(hb(ii),'facecolor',cols{ii},'barwidth',0.7,'EdgeColor','none','BaseValue',bv,'ShowBaseline','off');
%     end
    xlim([bins(1)-(incr/2) bins(end)+(incr/2)]);
    ylim([0 maxY]);
    set(gca,'TickDir','out','FontSize',7,'FontWeight','Normal','XTick',bins(1:2:end),'linewidth',0.25);
   
    if length(distD) > 1
        sigR = significanceTesting(distD);
    end

    if ~isempty(barGraph)
        barPos = barGraph{2};
        maxY = barGraph{4};
        ySpacing = barGraph{6};
        sigTestName = barGraph{8};
        sigColor = barGraph{10};
        baseValue = barGraph{12};
        pos = get(gca,'Position');
        axesPos = pos + [barPos(1) barPos(2) 0 0];
        axesPos(3:4) = [barPos(3) barPos(4)];
        hca = axes('Position',axesPos);hold on;
        hs = sigR.anova.multcompare.h; ps = sigR.anova.multcompare.p;
        plotBarsWithSigLines(sigR.means,sigR.sems,sigR.combs,[hs ps],'colors',cols,'sigColor',sigColor,'maxY',maxY,'ySpacing',ySpacing,'sigTestName',sigTestName,'sigLineWidth',0.25,'BaseValue',baseValue);
        xlim([0.4 0.6+length(sigR.means)]);
        set(gca,'TickDir','out','FontSize',5,'FontWeight','Normal');
    end

    if nargout == 4
        varargout{1} = sigR;
    end

%     axes(ha);
    if ~isempty(legs)
        putLegend(gca,legs,'colors',cols);
    end
    return;
end

if ismatrix(distD) && strcmp(do_mean,'Yes')
    for dd = 1:size(distD,2)
        allBars = [];
        for an = 1:size(distD,1)
            bd = distD{an,dd};
            if size(bd,1) > 1 || size(bd,2) > 1
                bd = bd(:);
            end
            [bar1 xs] = hist(bd,bins); bar1 = 100*bar1/sum(bar1);
%             [bar1 xs] = hist(bd,bins); bar1 = 100*bar1/length(bd);
            allBars = [allBars;bar1];
        end
        [mDist,semDist] = findMeanAndStandardError(cumsum(allBars,2));
%         [mDist,semDist] = findMeanAndStandardError(allBars);
        shadedErrorBar(bins,mDist,semDist,{'color',cols{dd},'linewidth',0.25},0.7);
    end
    ha = gca;
%     sigR = significanceTesting(distD);
    
%     if nargout == 4
%         varargout{1} = sigR;
%     end
    
    axes(ha);
    if ~isempty(legs)
        putLegend(gca,legs,'colors',cols);
    end
end