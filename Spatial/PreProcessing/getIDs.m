function ids = getIDs(raw,colA,selRows)
    Alphabets = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';
    col = strfind(Alphabets,colA);
    if ~exist('selRows','var')
        selRows = 1:size(raw,1);
    end
    for ii = 1:length(selRows)
        try
            ids(ii) = str2num((raw{selRows(ii),col}));
        catch
            ids(ii) = NaN;
        end
    end
end
