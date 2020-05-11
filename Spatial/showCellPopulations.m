function showCellPopulations


data = evalin('base','data');
mData = evalin('base','mData');
colors = mData.colors;
sigColor = mData.sigColor;
n = 0;
%%
selCells{1} = getCellPop('C_1_2_P');
selCells{2} = getCellPop('C_1_2_O');
selCells{3} = getCellPop('C_2');

% selCells{1} = getCellPop('C_2a_3_P');
% selCells{2} = getCellPop('C_2a_3_O');
% selCells{3} = getCellPop('C_3');

selCells{1} = getCellPop('C_3a_4_P');
selCells{2} = getCellPop('C_3a_4_O');
selCells{3} = getCellPop('C_4');

% selCells{1} = getCellPop('C_1');
% selCells{2} = getCellPop('C_2');
% selCells{3} = getCellPop('C_3');
% selCells{4} = getCellPop('C_4');

eii = 1;
cellInds = mData.cellsInds;
tei = evalin('base',sprintf('ei{%d}',eii));
% areCells = evalin('base',sprintf('ei{%d}.areCells',eii));
for ii = 1:length(selCells)
    thisList = cellInds(:,selCells{ii});
    inds = thisList(2,:)==eii;
    selCells{ii} = thisList(1,inds);
end

%%
figure(1000);clf
ha = axes;
showCells(ha,tei,selCells);

