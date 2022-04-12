function respB1 = distribution_of_cells_ocross(resp,bins)


binNs = unique(bins{1,1}); binNs = binNs(~isnan(binNs));
respB1 = [];
for rr = 1:size(resp,1)
    respB = [];
    for cc = 1:size(resp,2)
        tresp = resp(rr,cc);
        tplbin = bins{rr,cc};
        cell_resp = [];
        for bni = 1:length(binNs)
            cell_resp = [cell_resp cell_list_op({tplbin == binNs(bni)},tresp,'and')];
        end
        respB = [respB cell_resp];
    end
    respB1 = [respB1;respB];
end
