function figure_rasters(fn,allRs,ccs)

adata = evalin('base','dataT');
mData = evalin('base','mData');
colors = mData.colors;
sigColor = mData.sigColor;
mData.belt_length = adata{1}{1}{1}.belt_length;
selAnimals = 9;
n = 0;
%%


for ii = 1:2%length(data)
    distDi = [];
    SIs = [];rsi = [];pwsi = [];centersi = [];formulai = [];
    for jj = 1:length(selAnimals)
         [tempD cnsjj] = getVariableValues(adata{selAnimals(jj)},'placeCells5',ii);
         distDi = [distDi tempD];
         [tempS,~] = getVariableValues(adata{selAnimals(jj)},'SI',ii);
         SIs = [SIs tempS];
         [tempS,~] = getVariableValues(adata{selAnimals(jj)},'rs',ii);
         rsi = [rsi tempS];
         [tempS,~] = getVariableValues(adata{selAnimals(jj)},'centers',ii);
         centersi = [centersi tempS];
         [tempS,~] = getVariableValues(adata{selAnimals(jj)},'pws',ii);
         pwsi = [pwsi tempS];
         [tempS,~] = getVariableValues(adata{selAnimals(jj)},'formula',ii);
         formulai = [formulai tempS];
    end
    distD{ii} = distDi; SIsa{ii} = SIs; rss{ii} = rsi; centers{ii} = centersi; pws{ii} = pwsi;
    formula{ii} = formulai;
end
allPlaceCells = [];
for ii = 1:length(distD)
    selCells{ii} = find(distD{ii})
    allPlaceCells = [allPlaceCells selCells{ii}];
    oSelCells{ii} = setdiff(1:854,selCells{ii});
end
tpc = length(unique(allPlaceCells));
display(sprintf('Total number of place cells (unique) %d / %d = %.2f %%',tpc,length(distD{1}),100*tpc/length(distD{1})));
% selCells = selectCells(data,mData,'Common',[1 2]);
for ii = 1:length(distD)
    distDi = [];
    for jj = 1:length(selAnimals)
         [tempD cnsjj] = getVariableValues(adata{jj},'rasters',ii);
         distDi = cat(3,distDi,tempD);
    end
    data{ii}.rasters = distDi;
    data{ii}.dist = adata{1}{1}{1}.duration(1,:);
    data{ii}.SI = SIsa{ii};    data{ii}.rs = rss{ii}; data{ii}.centers = centers{ii};
    data{ii}.pws = pws{ii};     data{ii}.formula = formula{ii};
end
plotRastersMulti(data,1:size(data{1}.rasters,3),0,1,1)