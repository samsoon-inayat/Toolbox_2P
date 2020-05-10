function place_field_analysis_2_for_172110_Feb_02_2018_Recording_1

ei = evalin('base','ei{3}');
GF = gausswin(3);
b = ei.b;
signals = ei.tP.signals;

% S = ei.S;
% ccsi = find_place_cells(S,'adjRSquare_threshold',0.8);
ccsi = 1:length(ei.areCells);
ccs = ei.areCells(ccsi);

min_pc_width = 2;
max_pc_width = 120;
cm_per_bin = 150/100;

spSigAll = ei.deconv.spSigAll;
caSigAll = ei.deconv.caSigAll;
coi = [];
coiSI = [];
allSIs =[];
% load('temp.mat');
% labelCells(ei,coiSI,allSIs);
% OF_frames = OF_caSig(ei,ccsi(coi));
% for cn = 1:length(ccsi)

numberOfRows = 3;
numberOfCols = 3;
graphsOnOneFigure = numberOfRows * numberOfCols;
numberOfData = length(ccsi);
numberOfGroups = ceil(numberOfData/graphsOnOneFigure);
numberOfFullGroups = floor(numberOfData/graphsOnOneFigure);
indices = NaN(1,(numberOfGroups*graphsOnOneFigure));
indices(1:numberOfData) = 1:numberOfData;
groupIndices = reshape(indices,numberOfRows,numberOfCols,numberOfGroups);
for gg = 1:numberOfGroups
    groupIndices(:,:,gg) = groupIndices(:,:,gg)';
end

ff = makeFigureRowsCols(105,[NaN 0.5 11 4],'RowsCols',[numberOfRows numberOfCols],...
    'spaceRowsCols',[0.05 0.0009],'rightUpShifts',[0.03 0.03],'widthHeightAdjustment',...
    [-10 -55]);
ff_P = makeFigureRowsCols(106,[NaN 0.5 11 4],'RowsCols',[numberOfRows numberOfCols],...
    'spaceRowsCols',[0.05 0.0009],'rightUpShifts',[0.03 0.03],'widthHeightAdjustment',...
    [-10 -55]);

A_Trials = b.trials(1:20);
AP_Trials = b.trials(21:40);
trials = A_Trials;
onsets = b.air_puff_r(trials);
offsets = b.air_puff_f(trials);

photo_sensor = b.photo_sensor_f(b.photo_sensor_f>onsets(1) & b.photo_sensor_f<offsets(end));
diff_photo_sensor = diff(photo_sensor);
inds = find(diff_photo_sensor < 10000);
temp_photo_sensor = photo_sensor;
temp_photo_sensor(inds) = [];

a_onsets = temp_photo_sensor(1:(length(temp_photo_sensor)-1));
a_offsets = temp_photo_sensor(2:end);
% plotMarkers(ei,a_onsets,a_offsets,100);


trials = AP_Trials;
onsets = b.air_puff_r(trials);
offsets = b.air_puff_f(trials);

photo_sensor = b.photo_sensor_f(b.photo_sensor_f>onsets(1) & b.photo_sensor_f<offsets(end));
diff_photo_sensor = diff(photo_sensor);
inds = find(diff_photo_sensor < 100000);
temp_photo_sensor = photo_sensor;
temp_photo_sensor(inds) = [];

ap_onsets = temp_photo_sensor(1:(length(temp_photo_sensor)-1));
ap_offsets = temp_photo_sensor(2:end);
% plotMarkers(ei,ap_onsets,ap_offsets,100);

gg = 1;
while 1
    for rr = 1:numberOfRows
        for cc = 1:numberOfCols
            axes(ff.h_axes(rr,cc));
            cn = groupIndices(rr,cc,gg);
            if isnan(cn)
                cla(ff.h_axes(rr,cc));
                title('');
                continue;
            end
            tsp = spSigAll{ccsi(cn)}';
            caSig = signals(ccs(cn),:)';
            sCaSig = caSigAll{ccsi(cn)}';
            A = getDistRaster_1(b,caSig,tsp,a_onsets,a_offsets,0);
            plotDistRaster(ff.h_axes(rr,cc),A,[],num2str(cn));
            AP = getDistRaster_1(b,caSig,tsp,ap_onsets,ap_offsets,0);
            plotDistRaster(ff_P.h_axes(rr,cc),AP,[],num2str(cn));
        end
    end
    display('Press key');
    gg = keyboardInput(gg,[1 numberOfGroups],[1 6],'');
    if gg < 0
        break;
    end
end
