function data1 = combineData(data1,data2,fieldsToIgnore)
if isempty(data1)
    data1 = data2;
    return;
end

field_names = fields(data1);

if ~exist('fieldsToIgnore','var')
    fieldsToIgnore = [];
end

for ii = 1:length(field_names)
    thisField = field_names{ii};
    if sum(strcmp(fieldsToIgnore,thisField)) > 0
        continue;
    end
    try
        eval(sprintf('data1.%s = [data1.%s data2.%s];',thisField,thisField,thisField));
    catch
        eval(sprintf('data1.%s = [data1.%s;data2.%s];',thisField,thisField,thisField));
    end
end

