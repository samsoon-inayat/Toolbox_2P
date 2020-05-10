function rows = getSelectedCol(raw,colA,selRows)
    Alphabets = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';
    col = strfind(Alphabets,colA);
    if ~exist('selRows','var')
        selRows = 1:size(raw,1);
    end
    for ii = 1:length(selRows)
        if isnan(raw{selRows(ii),col})
%             if strcmp(colA,'A')
%                 rows(ii,1) = ' ';
%             else
                rows{ii,1} = ' ';
%             end
        else
%             if strcmp(colA,'A')
%                 rows(ii,1) = raw{selRows(ii),col};
%             else
                rows{ii,1} = raw{selRows(ii),col};
%             end
        end
    end
end
