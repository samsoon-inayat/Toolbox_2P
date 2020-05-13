function adjust_axis_label(hc,upos)

xlims = get(gca,'xlim');
pos = get(hc,'Position');
pos = pos + upos;
set(hc,'Position',pos);