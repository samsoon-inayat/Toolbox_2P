function [values,cns,cells] = getParamValues(varName,ei,planeNumbers,contextNumber,stimMarker,rasterType,selCells,maxDistTime)

if sum(strcmp(fields(ei{1}.b),'belt_length')) == 0
     tempB = get_mean_belt_length(ei,'protocol 15');
     ei{1}.b.belt_length = mean(tempB{1});
end

pos = strfind(varName,'.');
if ~isempty(pos)
    varNameParts{1,1} = varName(1:(pos(1)-1));
    if length(pos) == 1
        varNameParts{2,1} = varName((pos(1)+1):end);
    end
    for ii = 2:length(pos)
        varNameParts{ii,1} = varName((pos(ii-1)+1):(pos(ii)-1));
    end
    if length(pos) > 1
        varNameParts{ii+1,1} = varName((pos(ii)+1):end);
    end
else
    varNameParts{1,1} = varName;
end

values = [];
cns = [];
cells = logical([]);
for ee = 1:length(ei)
    tei = ei{ee};
    allplanes = tei.plane;
    if ~ischar(planeNumbers)
        if isnan(planeNumbers)
            planes = allplanes;
        else
            for ii = 1:length(planeNumbers)
                try
                    planes(ii) = allplanes(planeNumbers(ii));
                catch
                    continue;
                end
            end
        end
    else
        planes = allplanes;
    end
    if ~exist('planes','var')
        continue;
    end
    for pp = 1:length(planes)
        tplane = planes{pp};
        areCells = logical(tplane.tP.iscell(:,1));
        if ischar(selCells)
            if strcmp(selCells,'All')
                theSelectedCells = logical(ones(length(areCells),1));
            end
            if strcmp(selCells,'areCells')
                theSelectedCells = areCells;
            end
            if strcmp(selCells,'areCellsAll')
                theSelectedCells = logical(ones(sum(areCells),1));
            end
        else
            if size(selCells,2) == 1
                theSelectedCells = logical(selCells);
            else
                cellsInds = selCells(selCells(:,1) == pp,2);
                theSelectedCells = logical(zeros(size(areCells)));
                theSelectedCells(cellsInds) = 1;
            end
        end
        rasters = tplane.contexts(contextNumber).rasters;
        if strcmp(rasterType,'dist')
            cmdTxt = sprintf('data = rasters.%sD;',stimMarker);
        end
        if strcmp(rasterType,'time')
            cmdTxt = sprintf('data = rasters.%sT;',stimMarker);
        end
        eval(cmdTxt);
        tempV = []; tempC = []; 
        cells = [cells;areCells];
%         if strcmp(varNameParts{1},'dfbyfo')
%             cmdTxt = sprintf('data = rasters.%sT;',stimMarker);
%         end
        data.varName = varNameParts{1};
        if strcmp(varNameParts{1},'')
            values = populateDataProps(data,{rasterType,tei.b.belt_length,maxDistTime});
            cns = [];
%             values = combineData(values,tempV);
%             tempC = [ones(size(tempV.sp_rasters_nan_corrected,3),1)*ee ones(size(tempV.sp_rasters_nan_corrected,3),1)*pp find(theSelectedCells)];
%             cns = [cns;tempC];
        end
        if ~isempty(strfind(varNameParts{1},'data'))
            tempVZ = populateDataProps(data,{rasterType,tei.b.belt_length,maxDistTime});
            cmdTxt = sprintf('tempV = tempVZ.%s(theSelectedCells);',varNameParts{2});
            eval(cmdTxt);
             if isrow(tempV)
                tempV = tempV';
            end
            values = [values;tempV];
            tempC = [ones(size(tempV))*ee ones(size(tempV))*pp find(theSelectedCells)];
            cns = [cns;tempC];
        end
        if ~isempty(strfind(varNameParts{1},'cluster'))
            len_var_name = length(varNameParts{1});
            num_clus = str2double(varNameParts{1}(len_var_name));
            tempVZ = populateDataProps(data,{rasterType,tei.b.belt_length,maxDistTime,num_clus});
            cmdTxt = sprintf('tempV = tempVZ.%s & theSelectedCells;',varNameParts{1}(1:(len_var_name-1)));
            eval(cmdTxt);
             if isrow(tempV)
                tempV = tempV';
             end
            values = [values;tempV];
            tempC = [ones(size(tempV))*ee ones(size(tempV))*pp logical(ones(size(theSelectedCells)))];
            cns = [cns;tempC];
        end
        if ~isempty(strfind(varNameParts{1},'placeCells'))
            len_var_name = length(varNameParts{1});
            SI_th = str2double(varNameParts{1}(len_var_name));
            tempVZ = populateDataProps(data,{rasterType,tei.b.belt_length,maxDistTime,SI_th});
            cmdTxt = sprintf('tempV = tempVZ.%s(theSelectedCells);',varName);
%             cmdTxt = sprintf('tempV = tempVZ.%s(theSelectedCells);',varNameParts{1});
            eval(cmdTxt);
             if isrow(tempV)
                tempV = tempV';
             end
            tempV = logical(tempV);
            values = [logical(values);tempV];
            tempC = [ones(size(tempV))*ee ones(size(tempV))*pp find(theSelectedCells)];
            cns = [cns;tempC];
        end
        if strcmp(varNameParts{1},'rasters')
            tempVZ = populateDataProps(data,{rasterType,tei.b.belt_length,maxDistTime});
            cmdTxt = sprintf('tempV = tempVZ.rasters(:,:,theSelectedCells);',varName); eval(cmdTxt);
            values = cat(3,values,tempV);
            tempC = [ones(size(tempV,3),1)*ee ones(size(tempV,3),1)*pp find(theSelectedCells)];
            cns = [cns;tempC];
        end
        if strcmp(varNameParts{1},'fromFrames')
%             tempVZ = populateDataProps(data,{rasterType,tei.b.belt_length,maxDistTime});
            cmdTxt = sprintf('tempV = data.fromFrames.%s(:,:,theSelectedCells);',varNameParts{2}); eval(cmdTxt);
            if isempty(values)
                values = cat(3,values,tempV);
                tempC = [ones(size(tempV,3),1)*ee ones(size(tempV,3),1)*pp find(theSelectedCells)];
                cns = [cns;tempC];
            else
                svalues = size(values,2); stempV = size(tempV,2);
                if svalues > stempV
                    values = values(:,1:stempV,:);
                else
                    tempV = tempV(:,1:svalues,:);
                end
                values = cat(3,values,tempV);
                tempC = [ones(size(tempV,3),1)*ee ones(size(tempV,3),1)*pp find(theSelectedCells)];
                cns = [cns;tempC];
            end
        end
        if strcmp(varNameParts{1},'info_metrics') || strcmp(varNameParts{1},'fractal_dim') || strcmp(varNameParts{1},'place_field_properties')
            if strcmp(varNameParts{1},'place_field_properties')
                tempVZ = populateDataProps(data,{rasterType,tei.b.belt_length,maxDistTime});
                cmdTxt = sprintf('tempV = tempVZ.%s(theSelectedCells);',varName);
            else          
                cmdTxt = sprintf('tempV = data.%s(theSelectedCells);',varName);
            end
            eval(cmdTxt);
%             tempV = tempVZ.SI(theSelectedCells);
            if isrow(tempV)
                tempV = tempV';
            end
            values = [values;tempV];
            tempC = [ones(size(tempV))*ee ones(size(tempV))*pp find(theSelectedCells)];
            cns = [cns;tempC];
            
        end
        if strcmp(varNameParts{1},'gauss_fit_on_mean')
            if ~isempty(strfind(varName,'rs')) && isempty(strfind(varName,'rst'))
                tempVZ = populateDataProps(data,{rasterType,tei.b.belt_length,maxDistTime});
                cmdTxt = sprintf('tempV = tempVZ.%s(theSelectedCells);',varName); eval(cmdTxt);
                if isrow(tempV)
                    tempV = tempV';
                end
                values = [values;tempV];
                tempC = [ones(size(tempV))*ee ones(size(tempV))*pp find(theSelectedCells)];
                cns = [cns;tempC];
            end
            if ~isempty(strfind(varName,'rst'))
                tempVZ = populateDataProps(data,{rasterType,tei.b.belt_length,maxDistTime});
                cmdTxt = sprintf('tempV = tempVZ.%s(:,theSelectedCells);',varName); eval(cmdTxt);
                if isrow(tempV)
                    tempV = tempV';
                end
                values = [values tempV];
            end
            if ~isempty(strfind(varName,'coe'))
                cmdTxt = sprintf('tempV1 = data.%s.%s(theSelectedCells);','gauss_fit_on_mean',varNameParts{2}); eval(cmdTxt);
                for ii = 1:length(tempV1)
                    if ~isempty(tempV1(ii).rsquare)
                        tempV(ii,1) = tempV1(ii).rsquare;
                    else
                        tempV(ii,1) = NaN;
                    end
                end
                if isrow(tempV)
                    tempV = tempV';
                end
                values = [values;tempV];
                tempC = [ones(size(tempV))*ee ones(size(tempV))*pp find(theSelectedCells)];
                cns = [cns;tempC];
            end
        end
    end
end



function data = populateDataProps(data,inps)
% data.rasters = data.sp_rasters_nan_corrected;
binwidth = data.xs(2) - data.xs(1);
rasterType = inps{1};
belt_length = inps{2};
maxDistTime = inps{3};
if length(inps) == 4
    clus_SI_th = inps{4};
end
try
    mrfs = data.gauss_fit_on_mean;
catch
    if strcmp(rasterType,'dist')
        if maxDistTime == Inf
            binNum = size(data.sp_rasters_nan_corrected,2);
        else
            binNum = floor(maxDistTime(1)/binwidth);
        end
        data.rasters = data.sp_rasters_nan_corrected(:,1:binNum,:); % 48 corresponding to 141 cm
    else
        if maxDistTime == Inf
            binNum = size(data.sp_rasters_nan_corrected,2);
        else
            binNum = floor(maxDistTime(2)/binwidth);
        end
        data.rasters = data.sp_rasters_nan_corrected(:,1:binNum,:); % 21 corresponding to 6 secs
%         data.rasters = data.sp_rasters1;
    end
    return;
end
[rs,coeffs] = getMRFS_vals(mrfs);
data.gauss_fit_on_mean.rs = rs;
% [rst,coeffst] = getMRFS_vals_Trials(mrfs);
% data.gauss_fit_on_mean.rst = rst;
% data.gauss_fit_on_mean.coeffst = coeffst;
A = coeffs(:,1);
mu = coeffs(:,2);
sigma = coeffs(:,3);
as = A*exp(0.5);
bs = mu;
cs = sigma;
rs_th = 0.4;
PWs = 2.36*cs./sqrt(2)*binwidth;
data.place_field_properties.rs = rs; data.place_field_properties.pws = PWs';
data.place_field_properties.amp = as;
if strcmp(rasterType,'dist')
    data.place_field_properties.centers = (bs * binwidth)';
%     data.centers = (bs * binwidth)';
    if maxDistTime(1) == Inf
        binNum = size(data.sp_rasters_nan_corrected,2);
    else
        binNum = floor(maxDistTime(1)/binwidth);
    end
%     data.rasters = data.sp_rasters_nan_corrected(:,1:binNum,:); % 48 corresponding to 141 cm
    data.rasters = data.sp_rasters(:,1:binNum,:);
else
%     data.centers = (bs * binwidth)';
    data.place_field_properties.centers = (bs * binwidth)';
    if maxDistTime(2) == Inf
        binNum = size(data.sp_rasters_nan_corrected,2);
    else
        binNum = floor(maxDistTime(2)/binwidth);
    end
%     data.rasters = data.sp_rasters_nan_corrected(:,1:binNum,:); % 21 corresponding to 6 secs
    data.rasters = data.sp_rasters(:,1:binNum,:);
end
data.SI = data.info_metrics.ShannonMI_Zsh;

data.peaks = as';
center1 = size(data.rasters,1);
center2 = size(data.rasters,2);
if ~isempty(strfind(data.varName,'placeCells'))
%     cmdTxt = sprintf('data.placeCells%d = (logical((rs > rs_th) & (bs > center1 & bs < center2)'' & (data.SI > %d)))'';',clus_SI_th,clus_SI_th);
    cmdTxt = sprintf('data.placeCells%d = data.SI > %d;',clus_SI_th,clus_SI_th);
    eval(cmdTxt);
end


for ii = 1:size(data.sp_rasters,3)
    mRast(ii) = nanmean(nanmean(data.sp_rasters(:,:,ii)));
end
% data.propMatrix = [data.SI' data.fractal_dim.HaFD' data.fractal_dim.HiFD'];
% data.propMatrix = [data.SI' data.rs' data.fractal_dim.HaFD' data.fractal_dim.HiFD'];
data.propMatrix = [data.SI' data.fractal_dim.HiFD'];
% data.propMatrix = [data.SI' data.rs' data.fractal_dim.HaFD'];
% data.propMatrix = [data.SI' data.rs'];

rng('default');  % For reproducibility
if size(data.propMatrix,1)>1 & exist('clus_SI_th','var') & ~isempty(strfind(data.varName,'cluster'))
    rng(3,'twister');
%     eva = evalclusters(data.propMatrix,'kmeans','CalinskiHarabasz','KList',[1:25])
%     num_clus = eva.OptimalK;
    num_clus = 20;
    [clus_inds,clus_centers] = kmeans(data.propMatrix,num_clus);
%     data.cluster = clus_inds;
    for ii = 1:num_clus
        cluster_median_SI(ii) = median(data.SI(clus_inds == ii));
        cluster_mean_SI(ii) = mean(data.SI(clus_inds == ii));
        cluster_min_SI(ii) = min(data.SI(clus_inds == ii));
        cluster_min_HaFD(ii) = min(data.fractal_dim.HaFD(clus_inds == ii));
        cluster_min_HiFD(ii) = min(data.fractal_dim.HiFD(clus_inds == ii));
    end
%     [clus_inds2,clus_centers2] = kmeans([cluster_min_SI' cluster_min_HaFD' cluster_min_HiFD'],floor(num_clus/2));   
    clus_More = find(cluster_min_SI >= clus_SI_th);
%     clus_More = find(cluster_min_HaFD < 1);
%     clus_More = find(cluster_min_HiFD < 0.75);
    n_clus_inds = zeros(size(clus_inds));
    for ii = 1:length(clus_More)
        n_clus_inds = n_clus_inds | clus_inds == clus_More(ii);
    end
    data.cluster = n_clus_inds;% & rs' > 0.3;
%     data.cluster = data.fractal_dim.HaFD' < 0.7;
%     data.cluster = data.fractal_dim.HiFD' < 1.95 & data.SI' > clus_SI_th;
%     data.cluster = data.SI' > clus_SI_th;
end
n = 0;
% dur = data.duration_nan_corrected(:,1:93,:);
% rasters = data.rasters;
% parfor ii = 1:size(data.rasters,3)
%     thisRaster = rasters(:,:,ii);
%     tempCorr = corr(thisRaster);
%     tempCorr(isnan(tempCorr)) = 0;
%     corrs(:,:,ii) = tempCorr;
%     mtr = tril(tempCorr,-1);
% %     mcorr(ii) = mean(mtr(mtr ~= 0));
%     mcorr(ii) = BoxCountfracDim(mtr);
% end
% data.raster_trial_corr = corrs;
% data.mean_trial_corr = mcorr;


function data = combineData(data1,data2)
if isempty(data1)
    data = data2;
    return;
end

