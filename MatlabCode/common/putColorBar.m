function hc = putColorBar(hai,positionShift,minmax,fs,varargin)

p = inputParser;
addRequired(p,'hai',@isobject);
addRequired(p,'positionShift',@isnumeric);
addRequired(p,'minmax');
addRequired(p,'fs',@isnumeric);
addOptional(p,'location','eastoutside',@ischar);
addOptional(p,'text_margins',[0.3 0.125 0.3 0.125],@isnumeric);
parse(p,hai,positionShift,minmax,fs,varargin{:});

location = p.Results.location;
text_margins = p.Results.text_margins;

pos = get(hai,'Position');
ha = axes('Position',pos,'Visible','off');
axis off;
hc = colorbar('location',location);
set(hc,'EdgeColor','none');
changePosition(ha,positionShift);
% colormap parula;
if iscell(minmax)
    minmaxtextN{1} = minmax{1};
    minmaxtextP{1} = minmax{1};
    minmaxtextF{1} = minmax{1};
    minmaxtextN{2} = minmax{2};
    minmaxtextP{2} = minmax{2};
    minmaxtextF{2} = minmax{2};
    minmax = [0 1];
else
    minmaxtextN{1} = sprintf('-%d',-round(minmax(1)));
    minmaxtextP{1} = sprintf('%d',minmax(1));
    minmaxtextF{1} = sprintf('%.1f',minmax(1));
    
    minmaxtextN{2} = sprintf('-%d',-round(minmax(2)));
    minmaxtextP{2} = sprintf('%d',minmax(2));
    minmaxtextF{2} = sprintf('%.1f',minmax(2));
end
caxis(minmax);
xlims = xlim;
xincr = text_margins(1); yincr = text_margins(2);
xincr1 = text_margins(3); yincr1 = text_margins(4);
if strcmp(location,'northoutside')
    if minmax(1) > 10
        minmax(1) = round(minmax(1));
    end
    if minmax(2) > 10
        minmax(2) = round(minmax(2));
    end
    if abs(minmax(1)) > 10
        if minmax(1) < 0
            text(xlims(1)-xincr,1+yincr,minmaxtextN{1},'FontSize',fs);
        else
            text(xlims(1)-xincr,1+yincr,minmaxtextP{1},'FontSize',fs);
        end
    else
        text(xlims(1)-xincr,1+yincr,minmaxtextF{1},'FontSize',fs);
    end
    if abs(minmax(2)) > 10
        if minmax(2) < 0
            text(xlims(2)+xincr1,1+yincr1,minmaxtextN{2},'FontSize',fs);
        else
            text(xlims(2)+xincr1,1+yincr1,minmaxtextP{2},'FontSize',fs);
        end
    else
        text(xlims(2)+xincr1,1+yincr1,minmaxtextF{2},'FontSize',fs);
    end
    axis off; box off;
    set(hc,'Ticks',minmax,'TickLabels',{});
end
if strcmp(location,'eastoutside')
    if minmax(1) > 10
        minmax(1) = round(minmax(1));
    end
    if minmax(2) > 10
        minmax(2) = round(minmax(2));
    end
    if abs(minmax(1)) > 10
        if minmax(1) < 0
            text(xlims(2)+xincr,0-yincr,minmaxtextN{1},'FontSize',fs);
        else
            text(xlims(2)+xincr,0-yincr,minmaxtextP{1},'FontSize',fs);
        end
    else
        text(xlims(2)+xincr,0-yincr,minmaxtextF{1},'FontSize',fs);
    end
    if abs(minmax(2)) > 10
        if minmax(2) < 0
            text(xlims(2)+xincr1,1+yincr1,minmaxtextN{2},'FontSize',fs);
        else
            text(xlims(2)+xincr1,1+yincr1,minmaxtextP{2},'FontSize',fs);
        end
    else
        text(xlims(2)+xincr1,1+yincr1,minmaxtextF{2},'FontSize',fs);
    end
    axis off; box off;
    set(hc,'Ticks',minmax,'TickLabels',{});
end
minmax;


% try
%     set(hc,'Ticks',[minpcs maxpcs],'TickLabels',{sprintf('%.1f',(minpcs)),sprintf('%.1f',(maxpcs))},'FontSize',7);
% catch
% end