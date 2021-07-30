function [ha,hb,hca,varargout] = plotAverageDistributions (distD,varargin)
hb = NaN; hca = NaN; ha = NaN;
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
