function mdl = fit_linear(xs,ys)


statsetfitnlm = statset('fitnlm');
statsetfitnlm.MaxIter = 1000;
statsetfitnlm.TolFun = 1e-10;
% statsetfitnlm.Display = 'iter';
statsetfitnlm.TolX = statsetfitnlm.TolFun;
statsetfitnlm.UseParallel = 0;

mdl = fitnlm(xs,ys,'b(1)*x+b(2)',[0.5 0],'options',statsetfitnlm);
