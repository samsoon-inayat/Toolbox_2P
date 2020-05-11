function out = do_ancova(xa,ya,gr)

[handle,atab,ctab,stats1] = aoctool(xa,ya,gr,0.05,'Pred1','Response','Pred2','off','separate lines');
[mcc,mcm,mch,mcgnames] = multcompare(stats1,'Estimate','slope','CType','bonferroni','display','off');
out.ancova.anova_tab = atab;
out.ancova.stats = stats1;
out.ancova.multcompare.hs = mcc(:,6)<0.05;
out.ancova.multcompare.ps = mcc(:,6);
out.combs = mcc(:,1:2);
[mcc,mcm,mch,mcgnames] = multcompare(stats1,'Estimate','intercept','CType','bonferroni','display','off');
out.ancova.multcompare.hi = mcc(:,6)<0.05;
out.ancova.multcompare.pi = mcc(:,6);