function ff = show_remapping_corr_plots(mData,ff,remap_corrs,xs,minmaxC)

szRC = size(remap_corrs,1);
maskdisp = triu(ones(szRC,szRC),0);

FS = mData.axes_font_size;

for rr = 1:size(remap_corrs,1)
    for cc = 1:size(remap_corrs,2)
        minC(rr,cc) = min(remap_corrs{rr,cc}(:));
        maxC(rr,cc) = max(remap_corrs{rr,cc}(:));
    end
end

if ~isempty(minmaxC)
    m = minmaxC(1);
    M = minmaxC(2);
else
    m = min(minC(:));
    M = max(maxC(:));
end

for rr = 1:size(remap_corrs,1)
    for cc = 1:size(remap_corrs,2)
        if ~maskdisp(rr,cc)
            delete(ff.h_axes(rr,cc));
            continue;
        end
        axes(ff.h_axes(rr,cc));
        imagesc(remap_corrs{rr,cc},[m M]);
        box off;
        set(gca,'Ydir','Normal','linewidth',0.5,'FontSize',FS-1,'FontWeight','Bold');
        if rr == cc
            ts = xs.vals;
            set(gca,'YTick',xs.ticks,'YTickLabel',ts(xs.ticks));
            set(gca,'XTick',xs.ticks,'XTickLabel',ts(xs.ticks));
            xlabel(xs.label);
            ylabel(xs.label);
        else
            axis off
        end
        hc = putColorBar(ff.h_axes(rr,cc),[0.0 0.03 0 -0.05],[m M],5,'eastoutside',[0.07 0.08 0.1 0.13]);
    end
end