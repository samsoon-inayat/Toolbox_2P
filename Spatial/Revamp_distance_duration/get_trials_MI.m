function o = get_trials_MI(o)

n = 0;



nbins = 4;
nShuffle = 0;

Rs = o.Rs;
n = 0;

for rr = 1:size(Rs,1)
    for cc = 1:size(Rs,2)
        disp([rr cc]);
%         tic
        tRs = Rs{rr,cc};
        rast = tRs.sp_rasters1;
        numcells = size(rast,3);
        ntrials = size(rast,1);
        tMI = (NaN(size(rast,3),ntrials));
        for tni = 1:ntrials
            parfor cni = 1:numcells
                FR = squeeze(rast(tni,:,cni));
                tMI(cni,tni) = calc_MI(FR);
            end
        end
        tRs.MI_trials = tMI;
        Rs{rr,cc} = tRs;
%         toc
    end
end
o.Rs = Rs;


