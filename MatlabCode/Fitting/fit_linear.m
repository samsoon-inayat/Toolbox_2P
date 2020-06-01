function mdl = fit_linear(xs,ys)


statsetfitnlm = statset('fitnlm');
statsetfitnlm.MaxIter = 1000;
statsetfitnlm.TolFun = 1e-10;
% statsetfitnlm.Display = 'iter';
statsetfitnlm.TolX = statsetfitnlm.TolFun;
statsetfitnlm.UseParallel = 0;

mdl_fun = @(b,x)b(2)+b(1)*x;
mdl = fitnlm(xs,ys,mdl_fun,[0.5 0.5],'options',statsetfitnlm);
