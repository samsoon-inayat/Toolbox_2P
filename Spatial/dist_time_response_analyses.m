function dist_time_response_analyses

% contexts(1).name = 'Air Puff - Blank Belt - No Brake'; % There is no brake applied and animal runs in response to air puff for a fixed distance before air puff stops. No cues on the belt
% contexts(2).name = 'Air Puff - Cued Belt - No Brake'; % There is no brake applied and animal runs in response to air puff for a fixed distance before air puff stops. Cues on the belt
% contexts(3).name = 'Air Puff - Cued Belt - Brake'; % The air puff is locked to belt position using a brake. Cues on the belt
% contexts(4).name = 'Air Puff - Blank Belt - Brake'; % The air puff is locked to belt position using a brake. Blank belt
% contexts(5).name = 'Air Puff with brake - Blank Belt'; % Air trials only. The animal is stopped with brake. This is to see cellular responses to air puff
% contexts(6).name = 'All Air';
% contexts(7).name = 'Motion Intertrial';
% contexts(8).name = 'Motion Onsets';
% contexts(9).name = 'Motion Offsets';

thisCols_all = {[0 0 0],'b','r',[0 0.7 0.3],[0 0.7 0.3],'r',[0 0.7 0.3],[0 0 0],'b','m','r'};

aei = evalin('base','ei');
ei = aei([1 2 3 4 5 6]);
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
allCells = 1:sum(mData.cellsN);
% cells which respond to belt distance in context 1 and 2 --> landmark
% vector cells
sith = 10;
data_1 = a_ddata(1:5);
data_2 = b_ddata(1:5);
data = {data_1{1:2} data_2{1:4}};
ids = [1]; purePIs = selectCells(data,mData,ids,[1],[sith]); % these are cells which could be considered pure path integrators responding to reference as air onset
ids = [1 5]; notPurePIs = selectCells(data,mData,ids,[0 1],[sith sith]); % these are cells which are not pure path integrators necessarily
plotPopulationVectors_dual(106,data,data,{purePIs,notPurePIs},[1 5]);
plotPopulationSpatialInformation(107,data([1 2 5 6]),mData,allCells);
return;
spo = findPeaks_S(b_ddata{3},[],1);
% ids = [1 2 3 4]; lmvCells = selectCells(data,mData,ids,[1 0 1 1],[5 5 5 5]);
ids = [1]; purePIs = selectCells(data,mData,ids,[1],[5]); % these are cells which could be considered pure path integrators responding to reference as air onset
ids = [1 5]; notPurePIs = selectCells(data,mData,ids,[0 1],[5 5]); % these are cells which are not pure path integrators necessarily
ids = [2]; pureSensory = selectCells(data,mData,ids,[1],[7]);
ids = [1]; lmvCells = selectCells(moffaon_tdata(4),mData,ids,[1],[1]);
n = 0;
pfs = getPlaceCellPropsAuto(b_ddata(3),ids,lmvCells)
pfs = classifyCellsCNN(data,ids,lmvCells)

plotPopulationVectors(106,data,purePIs,1)
% ids = [8]; lmvCells = selectCells(mon_tdata,mData,ids,[1 1],[5]);
data = {b_ddata{3:4}};
dataT = {moffaon_tdata{3:4}};
plotPopulationVectors_dual(106,data,data,{purePIs,notPurePIs},[1 5]);
return;
plotPopulationVectors_dual(107,aon_tdata(6),aoff_tdata(6),lmvCells);
plotPopulationVectors_dual(108,mon_tdata(8),moff_tdata(9),lmvCells);


% plotPopulationSpatialInformation(107,data,mData,lmvCells);
return;
plotRastersMulti(data([1 2 5 6]),purePIs,0,0);

SI1 = data_bd{1}.SI(selCells{1});
SI2 = data_bd{2}.SI(selCells{1});
figure(106);clf;
scatter(SI1,SI2);
hold on;
SI1 = data_bd{1}.SI(selCells{3});
SI2 = data_bd{2}.SI(selCells{3});
scatter(SI1,SI2);

SI1 = a_ddata{6}.SI(selCells{1});
SI2 = a_tdata{6}.SI(selCells{1});
figure(106);clf;
scatter(SI1,SI2);




data = {a_tdata{1}};% ddata{2} tdata{2} ddata{3} tdata{3}};
data = {aon_tdata{4}};% ddata{2} tdata{2} ddata{3} tdata{3}};
plotPopulationVectors(106,data,selCells{1})
% 
% xlim([5 20]);
% ylim([5 20]);
% 

SI1 = mon_tdata{8}.SI(selCells{1});
SI2 = b_ddata{3}.SI(selCells{1});
figure(106);clf;scatter(SI1,SI2);hold on;
SI1 = mon_tdata{8}.SI(intersect(selCells{2},selCells{3}));
SI2 = b_ddata{3}.SI(intersect(selCells{2},selCells{3}));
scatter(SI1,SI2);

% plotPopulationSpatialInformation(107,[],data([2 4 6]),[],selCells{1});
% return;
% plotRastersMulti(data,selCells{1},0,0);

return
mData.cellsN = cellsN;
mData.sith = [10 10 10 10 5 5 5 5];
mData.dataList = dataList;
ind = 1;
for ii = 1:3:length(dataList)
    dataDef(ind).name = dataList{ii};
    dataDef(ind).trials = trials(dataList{ii+1}).def{2};
    dataDef(ind).functionName = dataList{ii+2};
    ind = ind + 1;
end

for ii = 1:length(dataDef)
    thisTrials = dataDef(ii).trials;
    cmdTxt = sprintf('func = @%s;',dataDef(ii).functionName);
    eval(cmdTxt);
    temp = func(tei,thisTrials);
    temp1 = find_mutual_information(tei,temp,owsi); temp.SI = temp1.zMI;
    temp.name = dataDef(ii).name;
    data{ii} = temp;
end

n=0;

ids = [3];
selCells{1} = selectCells(data,mData,ids,[1],[10]);
plotPopulationVectors(106,data(ids),selCells{1},1);
% plotCellPopulations(106,dataDef,data,mData)
plotPopulationSpatialInformation(107,dataDef,data([1 2 3 4]),mData,selCells{1});
selCells{1} = 1:sum(mData.cellsN);

pfs = getPlaceCellPropsGauss(brake_cued_belt,selCells{1});

% selCells_1 = segregateCells(selCells,cellsN,1);
% selCells_2 = segregateCells(selCells,cellsN,2);
% 
selCells_1 = segregateCells(pfs,cellsN,1);
showCells(107,tei{1},selCells_1);
showCells(108,tei{2},selCells_2);


% selCells{1} = find(brake_blank_belt.SI > 10 & no_brake_cued_belt.SI < 3 & motion_air_trials_onsets.SI < 3)
% selCells{1} = find(no_brake_cued_belt.SI > 7)




