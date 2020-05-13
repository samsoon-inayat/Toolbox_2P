% sample figure code

pp.numberOfRows = 2;
pp.spaceBetweenRows = [0.05];
pp.rowHeights = [0.15 0.45];

pp.numberOfCols = [1 1];
pp.colWidths = {[0.9]; [0.9]};
pp.spaceBetweenCols = {[0];[0]};

pp.leftOffset = 0.07;
pp.bottomOffset = 0.15;

pp = getPanelsProps(pp);


hf = makeFigureWindow(figNum,position,1);

pp = makeAxes(hf,pp,2,1,[0 0 0 0]);