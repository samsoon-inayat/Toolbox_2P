function makeSigTable(sigR,sigType,type,varargin)

p = inputParser;
default_colors = distinguishable_colors(20);
addRequired(p,'sigR');
addRequired(p,'sigType');
addRequired(p,'type');
addOptional(p,'position',[20 20 800 300],@isnumeric);
addOptional(p,'lines',0,@isnumeric);
parse(p,sigR,sigType,type,varargin{:});

position = p.Results.position;

sigD = sigR.combs;
if strcmp(sigType,'ks')
    sigD(:,3) = sigR.ks.p';
    sigD(:,4) = sigR.ks.h';
end

if strcmp(sigType,'annova')
    sigD(:,3) = sigR.annova.multcompare.p;
    sigD(:,4) = sigR.annova.multcompare.h;
end

if strcmp(sigType,'kruskalwallis')
    sigD(:,3) = sigR.kruskalwallis.multcompare.p;
    sigD(:,4) = sigR.kruskalwallis.multcompare.h;
end

if strcmp(type,'lines')
    diffCombs  = combs(:,2) - combs(:,1);
    maxdiff = max(combs(:,2) - combs(:,1));
    numberOfSigLines = length(find(sig(:,1)));
    
end

if strcmp(type,'table')
    pValTab = [];
    for ii = 1:size(sigD,1)
        this_p_val = sigD(ii,3);
        p_str = getNumberOfAsterisks(this_p_val);
        pValTab{sigD(ii,1),sigD(ii,2)} = sprintf('%s',p_str);
    end
    temp = sigD(:,[1 2]);
    gns = cellstr(num2str(sort(unique(temp(:)))));
    gns = gns';
    tbl = uitable(gcf,'Position',position);
    tbl.Data = pValTab;
    tbl.RowName = gns;
    tbl.ColumnName = gns;
    tbl.FontSize = 7;
    set(tbl,'ColumnWidth',{20 20 20 20});
end
