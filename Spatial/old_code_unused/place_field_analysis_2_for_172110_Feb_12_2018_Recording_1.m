function place_field_analysis_2_for_172110_Feb_02_2018_Recording_1

ei = evalin('base','ei{2}');
GF = gausswin(3);
b = ei.b;
signals = ei.tP.signals;

% S = ei.S;
% ccsi = find_place_cells(S,'adjRSquare_threshold',0.8);
ccsi = 1:length(ei.areCells);
ccs = ei.areCells(ccsi);

min_pc_width = 2;
max_pc_width = 120;
cm_per_bin = 150/50;

spSigAll = ei.deconv.spSigAll;
caSigAll = ei.deconv.caSigAll;
coi = [];
coiSI = [];
allSIs =[];
% load('temp.mat');
% labelCells(ei,coiSI,allSIs);
% OF_frames = OF_caSig(ei,ccsi(coi));
% for cn = 1:length(ccsi)

numberOfRows = 6;
numberOfCols = 2;
graphsOnOneFigure = numberOfRows * 1;
numberOfData = length(ccsi);
numberOfGroups = ceil(numberOfData/graphsOnOneFigure);
numberOfFullGroups = floor(numberOfData/graphsOnOneFigure);
indices = NaN(1,(numberOfGroups*graphsOnOneFigure));
indices(1:numberOfData) = 1:numberOfData;
groupIndices = reshape(indices,numberOfRows,1,numberOfGroups);
for gg = 1:numberOfGroups
    groupIndices(:,:,gg) = groupIndices(:,:,gg)';
end

ff = makeFigureRowsCols(105,[NaN 0.5 11 9],'RowsCols',[numberOfRows numberOfCols],...
    'spaceRowsCols',[0.05 0.045],'rightUpShifts',[0.03 0.03],'widthHeightAdjustment',...
    [-40 -55]);

A_Trials = b.trials(2:10);
AP_Trials = b.trials(13:23);
trials = A_Trials;
onsets = b.air_puff_r(trials);
offsets = b.air_puff_f(trials);
a_a_onsets = onsets;
a_a_offsets = offsets;

temp_on = b.air_puff_f([A_Trials]);
temp_off = b.air_puff_r([A_Trials]);
a_a_t_onsets = temp_on(1:(end-1));
a_a_t_offsets = temp_off(2:end);

temp_on = b.air_puff_f([AP_Trials]);
temp_off = b.air_puff_r([AP_Trials]);
ap_a_t_onsets = temp_on(1:(end-1));%[a_a_t_onsets;temp_on(1:(end-1))];
ap_a_t_offsets = temp_off(2:end);%[a_a_t_offsets;temp_off(2:end)];


photo_sensor = b.photo_sensor_f(b.photo_sensor_f>onsets(1) & b.photo_sensor_f<offsets(end));
diff_photo_sensor = diff(photo_sensor);
inds = find(diff_photo_sensor < 10000);
temp_photo_sensor = photo_sensor;
temp_photo_sensor(inds) = [];

a_p_onsets = temp_photo_sensor(1:(length(temp_photo_sensor)-1));
a_p_offsets = temp_photo_sensor(2:end);
% plotMarkers(ei,a_onsets,a_offsets,100);


trials = AP_Trials;
onsets = b.air_puff_r(trials);
offsets = b.air_puff_f(trials);
ap_a_onsets = onsets;
ap_a_offsets = offsets;

photo_sensor = b.photo_sensor_f(b.photo_sensor_f>onsets(1) & b.photo_sensor_f<offsets(end));
diff_photo_sensor = diff(photo_sensor);
inds = find(diff_photo_sensor < 100000);
temp_photo_sensor = photo_sensor;
temp_photo_sensor(inds) = [];

ap_p_onsets = temp_photo_sensor(1:(length(temp_photo_sensor)-1));
ap_p_offsets = temp_photo_sensor(2:end);
% plotMarkers(ei,ap_onsets,ap_offsets,100);

gg = 1;
while 1
    for rr = 1:numberOfRows
        cn = groupIndices(rr,1,gg);
        if isnan(cn)
            continue;
        end
        tsp = spSigAll{ccsi(cn)}';
        caSig = signals(ccs(cn),:)';
        sCaSig = caSigAll{ccsi(cn)}';
        A_P = getDistRaster_1(b,caSig,tsp,a_p_onsets,a_p_offsets,0); A_P_raster = A_P.distSigRaster./A_P.distDurRaster;
        A_P_SI = spatial_information(nanmean(A_P_raster),nanmean(A_P.distDurRaster));
        AP_P = getDistRaster_1(b,caSig,tsp,ap_p_onsets,ap_p_offsets,0); AP_P_raster = AP_P.distSigRaster./AP_P.distDurRaster;
        AP_P_SI = spatial_information(nanmean(AP_P_raster),nanmean(AP_P.distDurRaster));
        minD = min([A_P_raster(:);AP_P_raster(:)]);
        maxD = max([A_P_raster(:);AP_P_raster(:)]);
        for cc = 1:numberOfCols
            axes(ff.h_axes(rr,cc));
            if cc == 1
            	plotDistRaster(ff.h_axes(rr,cc),A_P,[minD maxD],sprintf('%s - %.2f',num2str(cn),A_P_SI));
            end
            if cc == 2
                plotDistRaster(ff.h_axes(rr,cc),AP_P,[minD maxD],sprintf('%s - %.2f',num2str(cn),AP_P_SI));
            end
%             if cc == 2
%                 Ar = getTimeRaster_simple(b,caSig,tsp,a_a_t_onsets,a_a_t_offsets,0);
%                 imagesc(Ar.sigRaster);colorbar;
%             end
%             if cc == 4
%                 Ar = getTimeRaster_simple(b,caSig,tsp,ap_a_t_onsets,ap_a_t_offsets,0);
%                 imagesc(Ar.sigRaster);colorbar;
%             end
        end
    end
    display('Press key');
    gg = keyboardInput(gg,[1 numberOfGroups],[1 6],'');
    if gg < 0
        break;
    end
end
