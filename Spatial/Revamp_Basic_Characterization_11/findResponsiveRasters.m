function out = findResponsiveRasters(rasters,Dur,ci,tRs,isCell)
p = NaN(size(rasters,3),1); p1 = p; p2 = p;
MIs = p;
if isempty(ci)
    MIs = tRs.info_metrics.ShannonMI_Zsh;
    out.zMIs = MIs;%(isCell);
    return;
end

if size(ci,2) == 3
group = [];
for ii = 1:size(ci,2)
    group = [group ii*ones(1,length(ci(1,ii):ci(2,ii)))];
end
% group = [ones(1,length(ci(1):ci(2))) (2*ones(1,(size(rasters,2)-ci)))];
p = NaN(size(rasters,3),1); p1 = p; p2 = p;
prm = p;
MIs = p;
trials = 1:size(rasters,1);
Dur = Dur(trials,:);
parfor ii = 1:size(rasters,3)
    thisRaster = rasters(trials,:,ii);
%     MIs(ii) = info_metrics_S_onlyMI_2(thisRaster,[],4,Dur,0);
%     [p(ii),atab,stats] = anova1(thisRaster,group,'nodisplay');
    [p(ii),atabk,statsk] = kruskalwallis(thisRaster,group,'nodisplay');
%     dataT = array2table(thisRaster);
%     within = array2table(group');
%     within.Var1 = categorical(within.Var1);
%     ra = repeatedMeasuresAnova(dataT,within,0.05);
%     prm(ii) = ra.ranova{3,5};
%     [p(ii),atab,stats] = ttest2(thisRaster,group,'nodisplay');
%     try
%     [mc,mm,mh,mgnames] = multcompare(stats,'CType','bonferroni','Display','off');
%     p1(ii) = mc(1,6);
%     p2(ii) = mc(3,6);
%     catch
%     end
%     if stats.means(1) < stats.means(2)
%         exin(ii) = 1;
%     end
%     if stats.means(1) > stats.means(2)
%         exin(ii) = -1;
%     end
%     if stats.means(1) == stats.means(2)
%         exin(ii) = 0;
%     end
%     thisRaster1 = thisRaster(:,ci(1,1):ci(2,1));
%     thisRaster2 = thisRaster(:,ci(1,2):ci(2,2));
%     thisRaster3 = thisRaster(:,ci(1,3):ci(2,3));
%     fd1(ii) = findHaFD(thisRaster1,1:size(thisRaster,1));
%     fd2(ii) = findHaFD(thisRaster2,1:size(thisRaster,1));
%     fd3(ii) = findHaFD(thisRaster3,1:size(thisRaster,1));
%     hifd1(ii) = findHiFD(thisRaster1,1:size(thisRaster,1));
%     hifd2(ii) = findHiFD(thisRaster2,1:size(thisRaster,1));
%     hifd3(ii) = findHiFD(thisRaster3,1:size(thisRaster,1));
%     if fd1(ii) < fd2(ii) & fd3(ii) < fd2(ii)
%         fdh(ii) = 1;
%     else
%         fdh(ii) = 0;
%     end
end
out.p = p;
% out.prm = prm;
% out.ps = [p1 p2];
out.h = p<0.05;
% out.MIs = MIs;
% out.exin = exin;
% out.fd1 = fd1;
% out.fd2 = fd2;
% out.fd3 = fd3;
% out.fdh = fdh;
% out.hifd1 = hifd1;
% out.hifd2 = hifd2;
% out.hifd3 = hifd3;
return;
end

if length(ci) == 1
group = [ones(1,ci) (2*ones(1,(size(rasters,2)-ci)))];
for ii = 1:size(rasters,3)
    thisRaster = rasters(:,:,ii);
    [p(ii),atab,stats] = anova1(thisRaster,group,'nodisplay');
    if stats.means(1) < stats.means(2)
        exin(ii) = 1;
    end
    if stats.means(1) > stats.means(2)
        exin(ii) = -1;
    end
    if stats.means(1) == stats.means(2)
        exin(ii) = 0;
    end
    thisRaster1 = thisRaster(:,1:ci);
    thisRaster2 = thisRaster(:,ci+1:end);
    fd1(ii) = findHaFD(thisRaster1,1:size(thisRaster,1));
    fd2(ii) = findHaFD(thisRaster2,1:size(thisRaster,1));
end
out.p = p;
out.h = p<0.05;
out.exin = exin;
out.fd1 = fd1;
out.fd2 = fd2;
return;
end