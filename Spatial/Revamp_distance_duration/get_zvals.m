function out = get_zvals(pop,ii)

out = pop;
for rr = 1:size(pop,1)
    for cc = 1:size(pop,2)
        zvals = pop{rr,cc};
        out{rr,cc} = zvals(:,ii);
    end
end

