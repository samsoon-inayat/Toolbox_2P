function selCellsO = segregateCells(selCells,cellsN,sel)

if sel == 1
    for ii = 1:length(selCells)
        selCellsO{ii} = intersect(selCells{ii},1:cellsN(1));
    end
elseif sel == 2
    for ii = 1:length(selCells)
        selCellsO{ii} = intersect(selCells{ii},(1:cellsN(2))+cellsN(1))-cellsN(1);
    end
end