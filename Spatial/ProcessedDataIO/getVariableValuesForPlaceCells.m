function [values,cns] = getVariableValuesForPlaceCells(data,varName,contextNumber,threshold)

[values,cns] = getVariableValues(data,varName,contextNumber);
[valuesSel,~] = getVariableValues(data,'sel',contextNumber);
[valuesSI,~] = getVariableValues(data,'SI',contextNumber);

inds = valuesSel & valuesSI > threshold;

values = values(inds);
cns = cns(:,inds);
