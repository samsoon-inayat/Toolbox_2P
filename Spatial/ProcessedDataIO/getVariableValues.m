function [values,cns] = getVariableValues(data,varName,contextNumber)

if ~strcmp(varName,'rasters') && ~strcmp(varName,'speed')
    values = [];
    cns = [];
    for ii = 1:length(data)
        tdata = data{ii};
        tdata = tdata{contextNumber};
        if strcmp(varName,'HaFD') || strcmp(varName,'HiFD')
            cmdTxt = sprintf('tempdata = tdata.fractal_dim.%s;',varName);
        else
            cmdTxt = sprintf('tempdata = tdata.%s;',varName);
        end
        eval(cmdTxt);
        cns = [cns [1:length(tempdata);(ones(size(tempdata))*ii)]];
        values = [values tempdata];
    end
    cns = [1:size(cns,2);cns];

%     inds = isnan(values);
%     values = values(~inds);
%     cns = cns(:,~inds);
return;
end
if strcmp(varName,'rasters')
    values = [];
    cns = [];
    for ii = 1:length(data)
        tdata = data{ii};
        tdata = tdata{contextNumber};
        cmdTxt = sprintf('tempdata = tdata.%s;',varName);
        eval(cmdTxt);
        cns = [cns [1:size(tempdata,3);(ones(1,size(tempdata,3))*ii)]];
        values = cat(3,values,tempdata);
    end
    cns = [1:size(cns,2);cns];
end

if strcmp(varName,'speed')
    values = [];
    cns = [];
    for ii = 1:length(data)
        tdata = data{ii};
        tdata = tdata{contextNumber};
        if ~isfield(tdata,'speedRasters')
            continue;
        end
        cmdTxt = sprintf('tempdata = tdata.%s;',varName);
        eval(cmdTxt);
        cns = [cns [1:size(tempdata,3);(ones(1,size(tempdata,3))*ii)]];
        values = cat(3,values,tempdata);
    end
    cns = [1:size(cns,2);cns];
end