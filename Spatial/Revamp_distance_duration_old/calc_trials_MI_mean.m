function o = calc_trials_MI_mean(o)

Rs = o.Rs;
n = 0;

for rr = 1:size(Rs,1)
    for cc = 1:size(Rs,2)
        disp([rr cc]);
        tRs = Rs{rr,cc};
        tRs.MI_trials_mean = nanmean(tRs.MI_trials);
        Rs{rr,cc} = tRs;
    end
end
o.Rs = Rs;


