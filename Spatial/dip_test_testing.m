xpdf = randn(1,10000);
xpdf = [xpdf randn(1,10000)+7];
[dip,xl,xu, ifault, gcm, lcm, mn, mj] = HartigansDipTest(xpdf);
[dip xl xu]

figure(100);clf;
subplot 211;
[vals,bins] = hist(xpdf,50);
plot(bins,vals);hold on;
plot([xl xl],[0 max(vals)],'r');
plot([xu xu],[0 max(vals)],'g');
subplot 212
vals = cumsum(vals/sum(vals));
plot(bins,vals); hold on;
% [vals1,bins1] = hist(gcm,50);
% plot(bins1,vals1);
title(sprintf('%d ifault, %.4f, %.4f, %.4f',ifault,dip,xl,xu));
plot([xl xl],[0 max(vals)],'r');
plot([xu xu],[0 max(vals)],'g');

