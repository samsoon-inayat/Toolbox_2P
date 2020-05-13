function set_axes_limits(ha,xlims,ylims)

if ~isempty(xlims)
    set(ha,'xlim',xlims);
end

if ~isempty(ylims)
    set(ha,'ylim',ylims);
end