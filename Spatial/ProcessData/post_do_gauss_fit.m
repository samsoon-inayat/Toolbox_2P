function c = post_do_gauss_fit (sig)

statsetfitnlm = statset('fitnlm');
statsetfitnlm.MaxIter = 1000;
statsetfitnlm.TolFun = 1e-10;
% statsetfitnlm.Display = 'iter';
statsetfitnlm.TolX = statsetfitnlm.TolFun;
statsetfitnlm.UseParallel = 1;
% statsetfitnlm.RobustWgtFun = 'welsch';
% fr = cell_act(mean_fr > 0.1,:);
bcs = 1:size(sig,2);
c = NaN(size(sig,1),4);
parfor ii = 1:size(sig,1)
    [~,~,c(ii,:)] = do_gauss_fit(bcs,sig(ii,:),statsetfitnlm,[1 0]);
end