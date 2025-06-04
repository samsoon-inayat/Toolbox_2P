
for ii = 1:5
    pd_rec = ei{ii};
    if ii == 1 || ii == 2
        trace_plot_all_dist(pd_rec,1); trace_plot_all_dist(pd_rec,2);
    else
        trace_plot_all_dist(pd_rec,1);
    end
end

