function display_colors(colors)

if ~exist('colors','var')
    colors = evalin('base','mData.shades.m');
end

ff = makeFigureRowsCols(107,[1 0.5 4 0.5],'RowsCols',[1 length(colors)],...
        'spaceRowsCols',[0 0.04],'rightUpShifts',[0.04 0.2],'widthHeightAdjustment',...
        [-50 -540]);
    set(gcf,'color','w');
    set(gcf,'Position',[5 3 10 1]);
    
for ii = 1:length(colors)
    ii
    ha = ff.h_axes(ii);
    set(ha,'Color',colors{ii});
%     axis off;
%     axis equal
    title(ii);
end
