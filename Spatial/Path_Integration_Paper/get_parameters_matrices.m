function out = get_parameters_matrices(aei,selAnimals,owr)

if ~exist('aei','var')
    aei = evalin('base','ei10');
    selAnimals = [1,2];
    owr = 1;
end

varNames = {'info_metrics.ShannonMI_Zsh','place_field_properties.amp','place_field_properties.pws','place_field_properties.centers','place_field_properties.rs'};
varNamesDH = {'zMIs','fFR','fwidths','fcenters','frs'};

if iscell(aei)

fileName = fullfile(pwd,'matFiles','parameters_matrics.mat');

if owr == 0
    out = load(fileName);
    return
end
% 

selCells = 'All';
planeNumbers = 'All';
maxDistTime = [142 15];


stimMarkers = {'air','air','belt','airI'};
rasterTypes = {'dist','time','dist','time'};
contextNumbers = [1 2 3 4];


for an = 1:length(selAnimals)
    tei = aei(selAnimals(an));
    for vi = 1:length(varNames)
        thisVarDH = varNamesDH{vi};
        cmdTxt = sprintf('%s_c = [];',thisVarDH);
        eval(cmdTxt);
    end
    for ci = 1:length(contextNumbers)
        contextNumber = contextNumbers(ci);
        for si = 1:length(stimMarkers)
            disp([an ci si]);
            stimMarker = stimMarkers{si};
            rasterType = rasterTypes{si};
            for vi = 1:length(varNames)
                thisVarDH = varNamesDH{vi};
                cmdTxt = sprintf('%s = [];',thisVarDH);
                eval(cmdTxt);
            end
            for vi = 1:length(varNames)
                thisVar = varNames{vi};
                thisVarDH = varNamesDH{vi};
                cmdTxt = sprintf('[%s cns areCells] = getParamValues(''%s'',tei,planeNumbers,contextNumber,stimMarker,rasterType,selCells,maxDistTime);',thisVarDH,thisVar);
                eval(cmdTxt);
                cmdTxt = sprintf('%s_c(ci,si,:) = %s;',thisVarDH,thisVarDH);
                eval(cmdTxt);
            end
        end
    end
    out.all_areCells{an} = areCells;
    out.all_cns{an} = cns;
    for vi = 1:length(varNames)
        thisVar = varNames{vi};
        thisVarDH = varNamesDH{vi};
        cmdTxt = sprintf('out.all_%s{an} = %s_c;',thisVarDH,thisVarDH);
        eval(cmdTxt);
    end
end

out.selCells = selCells;
out.planeNumbers = planeNumbers;
out.maxDistTime = maxDistTime;
out.stimMarkers = stimMarkers;
out.rasterTypes = rasterTypes;
out.contextNumbers = contextNumbers;

save(fileName,'-struct','out');
end

if isstruct(aei)
    paraMs = aei; clear aei;
    selC = selAnimals; clear selAnimals;
%     fieldNames = fields(selC);
    thresholdVars = {'zMI_threshold','','fwidth_limits','fcenter_limits','frs_threshold'};
    conditionNumbers = selC.conditionNumber;
    rasterTypes = selC.rasterType;
    for an = 1:length(paraMs.all_areCells)
        cellSel1s = logical(ones(size(paraMs.all_areCells{an})));
        cellSel0s = logical(zeros(size(paraMs.all_areCells{an})));
        cellSel = cellSel1s;
        if ~isnan(selC.areCells)
            if selC.areCells == 1
                cellSel = cellSel & paraMs.all_areCells{an};
            else
                cellSel = cellSel & ~paraMs.all_areCells{an};
            end
        end
        if ~isnan(selC.plane_number)
            cns = paraMs.all_cns{an};
            cellSel = cellSel & cns(:,2) == selC.plane_number;
        end
        for vi = 1:length(varNamesDH)
            if isempty(thresholdVars{vi})
                continue;
            end
            cmdTxt = sprintf('thisThreshold = selC.%s;',thresholdVars{vi});
            eval(cmdTxt);
            if isnan(thisThreshold)
                continue;
            end
            tempCR = cellSel0s;
            for ci = 1:length(conditionNumbers)
                CN = conditionNumbers(ci);
                tempC = cellSel0s;
                for ri = 1:length(rasterTypes)
                    RT = rasterTypes(ri);
                    cmdTxt = sprintf('tempR = squeeze(paraMs.all_%s{an}(CN,RT,:));',varNamesDH{vi});
                    eval(cmdTxt);
                    if length(thisThreshold) == 1
                        tempC = tempC | tempR > thisThreshold;
                    else
                        tempC = tempC | (tempR > thisThreshold(1) & tempR < thisThreshold(2));
                    end
                end
                tempCR = tempCR | tempC;
            end
            cellSel = cellSel & tempCR;
        end
%             if ~isnan(selC.zMI_threshold)
%                 temp = [];
%                 temp = squeeze(paraMs.all_zMIs{an}(selC.conditionNumber,selC.rasterType,:));
%                 cellSel = cellSel & temp > selC.zMI_threshold;
%             end
%             if ~isnan(selC.fcenter_limits)
%                 temp = [];
%                 temp = squeeze(paraMs.all_fcenters{an}(selC.conditionNumber,selC.rasterType,:));
%                 cellSel = cellSel & temp > selC.fcenter_limits(1) & temp < selC.fcenter_limits(2);
%             end
%             if ~isnan(selC.fwidth_limits)
%                 temp = [];
%                 temp = squeeze(paraMs.all_fwidths{an}(selC.conditionNumber,selC.rasterType,:));
%                 cellSel = cellSel & temp > selC.fwidth_limits(1) & temp < selC.fwidth_limits(2);
%             end
%             if ~isnan(selC.frs_threshold)
%                 temp = [];
%                 temp = squeeze(paraMs.all_frs{an}(selC.conditionNumber,selC.rasterType,:));
%                 cellSel = cellSel & temp > selC.frs_threshold;
%             end
        for vi = 1:length(varNamesDH)
            thisVarDH = varNamesDH{vi};
            cmdTxt = sprintf('out.all_%s{an} = paraMs.all_%s{an}(:,:,cellSel);',thisVarDH,thisVarDH);
            eval(cmdTxt);
        end
        out.cellSel{an} = cellSel;
    end
end

