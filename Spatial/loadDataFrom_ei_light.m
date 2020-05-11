function [data,mData] = loadDataFrom_ei_light(ei)

cei = combineRecordings(ei);
% cei.cellsN = size(ei{1}.signals,1);
% cei.cellsInds = 1:cei.cellsN;

mData.cellsN = cei.cellsN;
mData.cellsInds(1,:) = cei.cellsInds;
mData.cellsInds(2,:) = cei.cellsIndsRecs;
mData.sith = [5 5 5];
% mData.dataList = dataList;

selContexts = [3 4 5];
[a_ddata,a_tdata] = getDataContexts(cei,selContexts,'air');

mData.allCells = 1:sum(mData.cellsN);
data = a_ddata(1:3);
data = populate_mean_raster_fitting(data);

