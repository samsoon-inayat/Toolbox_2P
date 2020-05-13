function out = significanceTesting(data,varargin)
if isvector(data)
    %% ANOVA
    annovaVar = [];
    gr = [];
    for ii = 1:length(data)
        thisD = data{ii};
        thisDr = reshape(thisD,numel(thisD),1);
        annovaVar = [annovaVar;thisDr];
        gr = [gr;(ones(size(thisDr))*ii)];
        [mVals(ii) semVals(ii)] = findMeanAndStandardError(thisDr);
    end
    [p,tbl,stats] = anova1(annovaVar,gr,'off');%,names);
    out.annova.p = p;
    out.annova.tbl = tbl;
    out.annova.stats = stats;
    out.means = mVals;
    out.sems = semVals;
    % [p,tbl,stats] = kruskalwallis(y2wa,subst,'on');
    hf = figure(2001);
    [c,~,~,gnames] = multcompare(stats,'CType','hsd');
    out.annova.multcompare.c = c;
    out.annova.multcompare.gnames = gnames;
    close(hf);
    pdf = c(:,6);
    hdf = pdf<0.05;
    selGroups = 1:length(data);
    combs = nchoosek(1:length(selGroups),2);
    for ii = 1:size(combs,1)
        inds = ismember(c(:,[1 2]),[selGroups(combs(ii,1)) selGroups(combs(ii,2))],'rows');
        prT(ii,1) = c(inds,6);
    end
    hrT = prT<0.05;
    out.annova.multcompare.h = hrT;
    out.annova.multcompare.p = prT;
    out.combs = combs;


    %% Kruskal-Walis
    annovaVar = [];
    gr = [];
    for ii = 1:length(data)
        thisD = data{ii};
        thisDr = reshape(thisD,numel(thisD),1);
        annovaVar = [annovaVar;thisDr];
        gr = [gr;(ones(size(thisDr))*ii)];
    %     [mVals(ii) semVals(ii)] = findMeanAndStandardError(thisDr);
    end
    [p,tbl,stats] = kruskalwallis(annovaVar,gr,'off');%,names);
    out.kruskalwallis.p = p;
    out.kruskalwallis.tbl = tbl;
    out.kruskalwallis.stats = stats;
    % [p,tbl,stats] = kruskalwallis(y2wa,subst,'on');
    hf = figure(2001);
    [c,~,~,gnames] = multcompare(stats,'CType','hsd');
    out.kruskalwallis.multcompare.c = c;
    out.kruskalwallis.multcompare.gnames = gnames;
    close(hf);
    pdf = c(:,6);
    hdf = pdf<0.05;
    selGroups = 1:length(data);
    combs = nchoosek(1:length(selGroups),2);
    for ii = 1:size(combs,1)
        inds = ismember(c(:,[1 2]),[selGroups(combs(ii,1)) selGroups(combs(ii,2))],'rows');
        prT(ii,1) = c(inds,6);
    end
    hrT = prT<0.05;
    out.kruskalwallis.multcompare.h = hrT;
    out.kruskalwallis.multcompare.p = prT;
    return;
end

if ismatrix(data)
    if nargin > 1
        pointsToUse = varargin{1};
    end
    annovaVar = [];
    gr1 = [];
    gr2 = [];
    mVals = []; semVals = [];
    for ii = 1:size(data,1)
        for jj = 1:pointsToUse%size(data,2)
            thisD = data{ii,jj};
            thisDr = reshape(thisD,numel(thisD),1);
            if sum(isnan(thisDr))>0
                n = 0;
            end
            annovaVar = [annovaVar;thisDr];
            gr1 = [gr1;(ones(size(thisDr))*ii)];
            gr2 = [gr2;(ones(size(thisDr))*jj)];
            [mVals(ii,jj) semVals(ii,jj)] = findMeanAndStandardError(thisDr);
        end
    end
%     
%     groupType = ['1' '2' '3' '4']';
%     rmT = table(groupType,mVals(:,1),mVals(:,2),mVals(:,3),mVals(:,4),mVals(:,5),mVals(:,6)...
%     ,mVals(:,7),mVals(:,8),mVals(:,9),mVals(:,10),'VariableNames',{'Group','Center1','Center2','Center3',...
%     'Center4','Center5','Center6','Center7','Center8','Center9','Center10'});
%     Centers = [1:10]';
%     rm = fitrm(rmT,'Center1-Center10~Group','WithinDesign',Centers);
%     mauchlytbl = mauchly(rm)
%     ranovatbl = ranova(rm);
%     mcTbl = multcompare(rm,'Group','By','Time','ComparisonType','tukey-kramer')
%     mcTbl = multcompare(rm,'Group','ComparisonType','tukey-kramer')
%     mcTbl = multcompare(rm,'Time','ComparisonType','tukey-kramer')

    inds = isnan(annovaVar);
    [p,tbl,stats] = anovan(annovaVar,{gr1,gr2},'model','full','varnames',{'gr1','gr2'},'display','off');
    out.annova.p = p;
    out.annova.tbl = tbl;
    out.annova.stats = stats;
    out.means = mVals;
    out.sems = semVals;
    % [p,tbl,stats] = kruskalwallis(y2wa,subst,'on');
    hf = figure(2001);
    [c,~,~,gnames] = multcompare(stats,'CType','hsd','Dimension',[1 1]);
    out.annova.multcompare.c = c;
    out.annova.multcompare.gnames = gnames;
    close(hf);
    pdf = c(:,6);
    hdf = pdf<0.05;
    selGroups = 1:size(data,1);
    combs = nchoosek(1:length(selGroups),2);
    for ii = 1:size(combs,1)
        inds = ismember(c(:,[1 2]),[selGroups(combs(ii,1)) selGroups(combs(ii,2))],'rows');
        prT(ii,1) = c(inds,6);
    end
    hrT = prT<0.05;
    out.annova.multcompare.h = hrT;
    out.annova.multcompare.p = prT;
    out.combs = combs;
    n = 0;
end