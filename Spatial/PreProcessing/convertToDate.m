function dateC = convertToDate(inputDates)


for ii = 1:length(inputDates)
    thisDate = inputDates(ii);
    pos = strfind(thisDate,'/');
    if isempty(pos{1})
        dateC(ii,1) = datetime(thisDate,'Format','yyyy-MM-dd');
    else
        dateC(ii,1) = datetime(thisDate,'Format','MM/dd/yyyy');
    end
end