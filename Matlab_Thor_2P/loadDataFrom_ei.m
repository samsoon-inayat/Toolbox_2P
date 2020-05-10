function [data,mData] = loadDataFrom_ei(ei)

cei = combineRecordings(ei);
% cei.cellsN = size(ei{1}.signals,1);
% cei.cellsInds = 1:cei.cellsN;

mData.cellsN = cei.cellsN;
mData.cellsInds(1,:) = cei.cellsInds;
mData.cellsInds(2,:) = cei.cellsIndsRecs;
mData.sith = [10 10 10 10 5 5 5 5];
% mData.dataList = dataList;

selContexts = [1 2 3 4 5 6 7 8 9];
[a_ddata,a_tdata] = getDataContexts(cei,selContexts,'air');
[ai_ddata,ai_tdata] = getDataContexts(cei,selContexts,'airI');
[moffaon_ddata,moffaon_tdata] = getDataContexts(cei,selContexts,'motionOffsetAirOnset');
[b_ddata,b_tdata] = getDataContexts(cei,selContexts,'belt');
[moo_ddata,moo_tdata] = getDataContexts(cei,selContexts,'motionOnsetsOffsets');
[m_ddata,m_tdata] = getDataContexts(cei,selContexts,'motionI');
[mon_ddata,mon_tdata] = getDataContexts(cei,selContexts,'motionOnsets22');
[moff_ddata,moff_tdata] = getDataContexts(cei,selContexts,'motionOffsets22');
[aon_ddata,aon_tdata] = getDataContexts(cei,selContexts,'airOnsets22');
[aoff_ddata,aoff_tdata] = getDataContexts(cei,selContexts,'airOffsets22');

mData.allCells = 1:sum(mData.cellsN);
data_1 = a_ddata(1:5);
data_2 = b_ddata(1:5);
data = {data_1{1:2} data_2{1:4}};
data = populate_mean_raster_fitting(data);
data = data([1 2 5 6]);
