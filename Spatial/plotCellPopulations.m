function plotCellPopulations(fn,dataDef,data,mData)

n = 0;
totalCells = sum(mData.cellsN);
sith = mData.sith;
for ii = 1:length(data)
    cellPop = zeros(length(data),totalCells);
    selCells{1} = find(data{ii}.SI > sith(ii));
    ids = 1:length(data);
    ids(ii) = [];
    for jj = 1:length(ids)
       selCells{jj+1} = find(data{ii}.SI > 5 & data{ids(jj)}.SI > sith(ids(jj)));
    end
    for jj = 1:length(selCells)
        cellPop(jj,selCells{jj}) = 1*jj;
    end
    figure(fn);clf;
    imagesc(cellPop);colorbar;
    set(gca,'Ydir','Normal');
    n = 0;
end
