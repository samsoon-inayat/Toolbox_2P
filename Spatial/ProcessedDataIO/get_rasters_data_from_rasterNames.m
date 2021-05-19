function rasters = get_rasters_data_from_rasterNames(ei,rasterNames)
if ~exist('ei','var')
    ei = evalin('base','ei');
    rasterNames = {'airD','beltD'};
end

rasters = cell(length(ei),length(rasterNames));
for rr = 1:length(rasterNames)
    for ii = 1:length(ei)
        rasters{ii,rr} = get_data(ei{ii},rasterNames{rr});
    end
end


function rasters = get_data(ei,rasterName)
rasters = {};
contextNs = getContexts(ei,rasterName);
if isempty(contextNs)
    return;
else
    for ii = 1:length(contextNs)
        rasters{ii} = get_context_rasters(ei,contextNs(ii),rasterName)
        rasters{ii}.contextInfo = {contextNs(ii),ei.plane{1}.contexts(contextNs(ii)).name};
    end
end

function contextNs = getContexts(ei,rasterName)
pp = 1;
contexts = ei.plane{1}.contexts;
contextNs = [];
for ii = 1:length(contexts)
    these_rasters = contexts(ii).rasters;
    fields = fieldnames(these_rasters);
    if any(strcmp(fields,rasterName))
        contextNs = [contextNs ii];
    end
end


function rasters = get_context_rasters(ei,contextN,rasterName)
nplanes = length(ei.plane);
pp = 1;
thisContext = ei.plane{pp}.contexts(contextN);
disp(sprintf('%s - plane-%d',thisContext.name,pp));
cmdTxt = sprintf('tempR = thisContext.rasters.%s;',rasterName);
eval(cmdTxt);
iscell1 = ei.plane{pp}.tP.iscell(:,1);
if nplanes == 1
    rasters = tempR;
    rasters.iscell = iscell1;
else
    pp = 2;
    thisContext = ei.plane{pp}.contexts(selContexts(ii));
    disp(sprintf('%s - plane-%d',thisContext.name,pp));
    cmdTxt = sprintf('tempR1 = thisContext.rasters.%s;',rasterNames{ii});
    eval(cmdTxt);
    rasters = combine_planes_data(tempR,tempR1);
    iscell2 = ei.plane{pp}.tP.iscell(:,1);
    rasters.iscell = cat(1,iscell1,iscell2);
end


function tempRC = combine_planes_data(tempR,tempR1)
struct_fields_to_combine = {'fromFrames','info_metrics','gauss_fit_on_mean','fractal_dim'};
fields = fieldnames(tempR);
cellmatflag = 1;
for ff = 1:length(fields)
    thisFieldTxt = fields{ff};
    cmdTxt = sprintf('thisField = tempR.%s;',fields{ff});
    eval(cmdTxt);
    cmdTxt = sprintf('thisField1 = tempR1.%s;',fields{ff});
    eval(cmdTxt);
    if isstruct(thisField)
        if strcmp(fields{ff},'wholeData')
            thisFieldC = combine_whole_data(thisField,thisField1);
        end
        if strcmp(fields{ff},'fromFrames')
            thisFieldC = combine_planes_data(thisField,thisField1);
        end
        if strcmp(fields{ff},'info_metrics')
            thisFieldC = combine_info_metrics(thisField,thisField1);
        end
        if strcmp(thisFieldTxt,'gauss_fit_on_mean')
            thisFieldC = combine_gauss_fit(thisField,thisField1);
        end
        if strcmp(thisFieldTxt,'fractal_dim')
            thisFieldC = combine_fractal_dim(thisField,thisField1);
        end
        n = 0;
    else
        if size(thisField,3) == 1
            thisFieldC = thisField;
        else
            try
                thisFieldC = cat(3,thisField,thisField1);
            catch
                szf1 = size(thisField,2); szf2 = size(thisField1,2);
                df = szf1 - szf2;
                if df == 0
                    error;
                end
                if df > 0
                    thisFieldC = cat(3,thisField(:,1:szf2,:),thisField1);
                end
                if df < 0
                    thisFieldC = cat(3,thisField,thisField1(:,1:szf1,:));
                end
            end
        end
        if cellmatflag
            cell_numbers = [[1:size(thisFieldC,3)]' [ones(size(thisField,3),1);(2*ones(size(thisField1,3),1))] [[1:size(thisField,3)]';[1:size(thisField1,3)]']];
            cellmatflag = 0;
        end
    end
    cmdTxt = sprintf('tempRC.%s = thisFieldC;',fields{ff});
    eval(cmdTxt);
end
tempRC.cell_numbers = cell_numbers;

function ifC = combine_info_metrics(if1,if2)
% fields = fieldnames(if1);
% for ff = 1:length(fields)
%     cmdTxt = sprintf('ifC.%s = [];',fields{ff}); eval(cmdTxt);
% end
ifC.ShannonMI = cat(2,if1.ShannonMI,if2.ShannonMI);
ifC.ShannonMI_p = cat(2,if1.ShannonMI_p,if2.ShannonMI_p);
ifC.ShannonMI_Zsh = cat(2,if1.ShannonMI_Zsh,if2.ShannonMI_Zsh);
ifC.ShannonMI_shuffle = cat(1,if1.ShannonMI_shuffle,if2.ShannonMI_shuffle);


function thisFieldC = combine_gauss_fit(thisField,thisField1)
thisFieldC.coefficients_Rs_mean = cat(1,thisField.coefficients_Rs_mean,thisField1.coefficients_Rs_mean);
thisFieldC.worked = cat(1,thisField.worked,thisField1.worked);
thisFieldC.gauss1Formula = thisField.gauss1Formula;

function thisFieldC = combine_fractal_dim(thisField,thisField1)
thisFieldC.HaFD = cat(2,thisField.HaFD,thisField1.HaFD);
thisFieldC.HiFD = cat(2,thisField.HiFD,thisField1.HiFD);


function thisFieldC = combine_whole_data(thisField,thisField1)
fields = fieldnames(thisField);
for ii = 1:length(fields)
    fieldTxt = fields{ii};
    if strcmp(fieldTxt,'sp_rasters')
        thisFieldC.sp_rasters = cat(2,thisField.sp_rasters,thisField1.sp_rasters);
    elseif strcmp(fieldTxt,'cell_history')
        thisFieldC.cell_history = cat(1,thisField.cell_history,thisField1.cell_history);
    else
        cmdTxt = sprintf('thisFieldC.%s = thisField.%s;',fields{ii},fields{ii});
        eval(cmdTxt);
    end
end


