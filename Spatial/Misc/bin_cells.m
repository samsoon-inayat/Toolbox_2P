function [popS,perc1] = bin_cells(PL,bins,pop)
perc = [];
n = 0;
for rr = 1:size(PL,1)
    ind = 1;
    tp = [];
    for cc = 1:size(PL,2)
        bs = PL{rr,cc};
        [N,E,Bi] = histcounts(bs,bins);
        for nn = 1:length(N)
            popS{rr,ind} = logical((Bi == nn)) & pop{rr,cc};
            perc1(rr,ind) = 100*sum(popS{rr,ind})/length(bs);
            ind = ind + 1;
        end
        tp = [tp,100*N/length(bs)];
    end
    perc = [perc;tp];
end
