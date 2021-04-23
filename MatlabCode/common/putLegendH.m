function putLegend(ah,legs,varargin)
% function putLegend(ah,legs,specs,varargin)

p = inputParser;
default_colors = distinguishable_colors(20);
default_ySpacingFactor = 10;
addRequired(p,'ah',@ishandle);
addRequired(p,'legs',@iscell);
addOptional(p,'colors',default_colors,@iscell);
addOptional(p,'sigR',NaN);
addOptional(p,'sigType',NaN);
addOptional(p,'sigColor',NaN);
addOptional(p,'lineWidth',0.5);


parse(p,ah,legs,varargin{:});

cols = p.Results.colors;
lineWidth = p.Results.lineWidth;
ah = p.Results.ah;
temp = p.Results.sigR;
if iscell(temp)
    sigR = temp{1};
    sigType = temp{2};
    sigColor = temp{3};
    sigFontSize = temp{4};
else
    sigR = NaN;
    sigFontSize = 5;
end

axes(ah);
specs = legs{end};
legs(end) = [];
lenLine = specs(2);
x1 = specs(1); x2 = x1+lenLine; y1 = specs(3); gap = specs(4);
legendFontSize = sigFontSize;
for ii = 1:length(legs)
    if isempty(legs{ii})
        continue;
    end
    plot([x1 x2],[y1 y1],'color',cols{ii},'linewidth',lineWidth);
    ht = text(x2+(lenLine/2),y1,sprintf('%s',legs{ii}),'Color',cols{ii},'FontSize',legendFontSize);
    te = get(ht,'Extent');
    x1 = te(1)+te(3)+gap;
    x2 = x1 + lenLine;
end

