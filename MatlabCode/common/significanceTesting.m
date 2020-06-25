function out = significanceTesting(data,varargin)

% p = inputParser;
% addRequired(p,'data');
% addOptional(p,'maxY',default_maxY,@isnumeric);
% 
% parse(p,means,sems,combs,sig,varargin{:});

%%
if isvector(data)
    anovaVar = [];
    gr = [];
    for ii = 1:length(data)
        thisD = data{ii};
        thisDr = reshape(thisD,numel(thisD),1);
        anovaVar = [anovaVar;thisDr];
        gr = [gr;(ones(size(thisDr))*ii)];
        [mVals(ii) semVals(ii)] = findMeanAndStandardError(thisDr);
    end
    out = doAnovaAndKruskalWallis(anovaVar,gr);
    out.means = mVals;
    out.sems = semVals;
    
    selGroups = 1:length(data);
    combs = nchoosek(1:length(selGroups),2);
    for ii = 1:size(combs,1)
        data1 = data{selGroups(combs(ii,1))};
        data2 = data{selGroups(combs(ii,2))};
        try
            [htt2(ii),ptt2(ii),citt2(ii,:),statstt2(ii)] = ttest2(reshape(data1,1,numel(data1)),reshape(data2,1,numel(data2)));
        catch
        end
        try
            [htt(ii),ptt(ii),citt(ii,:),statstt(ii)] = ttest(reshape(data1,1,numel(data1)),reshape(data2,1,numel(data2)));
        catch
        end
        try
            [hks(ii),pks(ii),ks2stat(ii)] = kstest2(reshape(data1,1,numel(data1)),reshape(data2,1,numel(data2)));
        catch
        end
    end
    out.combs = combs;
    if exist('htt2','var')
        out.ttest2.h = htt2; out.ttest2.p = ptt2; out.ttest2.ci = citt2; out.ttest2.stats = statstt2;
    end
    if exist('htt','var')
        out.ttest.h = htt; out.ttest.p = ptt; out.ttest.ci = citt; out.ttest.stats = statstt;
    end
    if exist('hks','var')
        out.kstest.h = hks; out.kstest.p = pks; out.kstest.ks2stat = ks2stat;
    end
    return;
end
%%
if ismatrix(data)
    if ~iscell(data)
        for ii = 1:size(data,2)
            mdatac{ii} = data(:,ii);
        end
        out = significanceTesting(mdatac);
        return;
    end
    if iscell(data)
        if isvector(data{1})
            for ii = 1:size(data,1) % go through each row
                try
                    rowSigR{ii} = significanceTesting(data(ii,:));
                catch
                    rowSigR{ii} = [];
                    disp('Row Significance test didn''t work');
                end
                for jj = 1:size(data,2)
                    mdata(ii,jj) = nanmean(data{ii,jj});
                end
            end
            for ii = 1:size(data,2)
                mdatac{ii} = mdata(:,ii);
            end
            out = significanceTesting(mdatac);
            out.rowSigR = rowSigR;
            
            dataRMA = mdata; numCols = size(dataRMA,2);
            cmdTxt = sprintf('dataT = table(');
            for ii = 1:(size(dataRMA,2)-1)
                cmdTxt = sprintf('%sdataRMA(:,%d),',cmdTxt,ii);
            end
            cmdTxt = sprintf('%sdataRMA(:,size(dataRMA,2)));',cmdTxt);
            eval(cmdTxt);
            varNames = dataT.Properties.VariableNames;
            cmdTxt = sprintf('rm = fitrm(dataT,''');
            for ii = 1:(length(varNames)-1)
                cmdTxt = sprintf('%s%s,',cmdTxt,varNames{ii});
            end
            cmdTxt = sprintf('%s%s~1'');',cmdTxt,varNames{length(varNames)});
            eval(cmdTxt);
            out.ranova.rm = rm;
            out.ranova.table = ranova(rm);
            out.ranova.mauchlytbl = mauchly(rm);
            mcTI = find_sig_mctbl(multcompare(rm,'Time','ComparisonType','bonferroni'),5);
            [combs,h,p] = populate_multcomp_h_p(dataRMA,rm.WithinDesign,mcTI,[]);
            out.ranova.multcompare.tbl = mcTI;
            out.ranova.multcompare.combs = combs;
            out.ranova.multcompare.h = h;
            out.ranova.multcompare.p = p;
            out.ranova.ds = descriptiveStatistics(dataRMA);
            return;
        end
    end    
%     if nargin > 1
%         pointsToUse = varargin{1};
%     else
%         pointsToUse = size(data,2);
%     end
%     anovaVar = [];
%     gr1 = [];
%     gr2 = [];
%     mVals = []; semVals = [];
%     for ii = 1:size(data,1)
%         for jj = 1:pointsToUse%size(data,2)
%             thisD = data{ii,jj};
%             thisDr = reshape(thisD,numel(thisD),1);
%             if sum(isnan(thisDr))>0
%                 n = 0;
%             end
% %             anovaVar = [anovaVar;thisDr];
%             gr1 = [gr1;(ones(size(thisDr))*ii)];
%             gr2 = [gr2;(ones(size(thisDr))*jj)];
%             [mVals(ii,jj) semVals(ii,jj)] = findMeanAndStandardError(thisDr);
%         end
%     end
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
    anovaVar = mVals';
    [p,tbl,stats] = anovan(anovaVar,1:size(anovaVar,2),'model','full','varnames',{'gr1','gr2'},'display','off');
    out.anova.p = p;
    out.anova.tbl = tbl;
    out.anova.stats = stats;
    out.means = mVals;
    out.sems = semVals;
    % [p,tbl,stats] = kruskalwallis(y2wa,subst,'on');
    hf = figure(2001);
    [c,~,~,gnames] = multcompare(stats,'CType','bonferroni','Dimension',[1 1]);
    out.anova.multcompare.c = c;
    out.anova.multcompare.gnames = gnames;
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
    out.anova.multcompare.h = hrT;
    out.anova.multcompare.p = prT;
    out.combs = combs;
    out.gr1 = gr1;
    out.gr2 = gr2;
    
    anovaVar = []; gr1 = []; gr2 = [];
    for ii = 1:size(data,1)
        for jj = 1:size(data,2)
            thisD = data{ii,jj};
            thisDr = reshape(thisD,numel(thisD),1);
            if sum(isnan(thisDr))>0
                n = 0;
            end
            anovaVar = [anovaVar;thisDr];
            gr1 = [gr1;(ones(size(thisDr))*ii)];
            gr2 = [gr2;(ones(size(thisDr))*jj)];
        end
    end
    
    
    [handle,atab,ctab,stats1] = aoctool(gr2,anovaVar,gr1,0.05,'Pred1','Var','Pred2','off','separate lines');
    [mcc,mcm,mch,mcgnames] = multcompare(stats1,'Estimate','slope','CType','bonferroni','display','off');
    out.ancova.anova_tab = atab;
    out.ancova.stats = stats1;
    out.ancova.multcompare.h = mcc(:,6)<0.05;
    out.ancova.multcompare.p = mcc(:,6);
    n = 0;
end


function out = doAnovaAndKruskalWallis(anovaVar,gr)
[p,tbl,stats] = anova1(anovaVar,gr,'off');%,names);
out.anova.p = p;
out.anova.tbl = tbl;
out.anova.stats = stats;
hf = figure(2001);
[c,~,~,gnames] = multcompare(stats,'CType','bonferroni');
out.anova.multcompare.c = c;
out.anova.multcompare.gnames = gnames;
close(hf);
ps = c(:,6);
hs = ps<0.05;
out.anova.multcompare.combs = c(:,[1 2]);
out.anova.multcompare.h = hs;
out.anova.multcompare.p = ps;

[p,tbl,stats] = kruskalwallis(anovaVar,gr,'off');%,names);
out.kruskalwallis.p = p;
out.kruskalwallis.tbl = tbl;
out.kruskalwallis.stats = stats;
% [p,tbl,stats] = kruskalwallis(y2wa,subst,'on');
hf = figure(2001);
[c,~,~,gnames] = multcompare(stats,'CType','bonferroni');
out.kruskalwallis.multcompare.c = c;
out.kruskalwallis.multcompare.gnames = gnames;
close(hf);
ps = c(:,6);
hs = ps<0.05;
out.anova.multcompare.combs = c(:,[1 2]);
out.kruskalwallis.multcompare.h = hs;
out.kruskalwallis.multcompare.p = ps;