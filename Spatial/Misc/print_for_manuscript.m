function print_for_manuscript(ra,row)

if exist('row','var')
    if isnumeric(row)
        ind = row;
    else
        ind = row_identifier(ra.ranova.Row,row);
    end
    DF1 = ra.ranova.DF(ind); DF2 = ra.ranova.DF(ind+1);
    F = ra.ranova.F(ind); 
    try
        p = ra.ranova.pValue_sel(ind); 
    catch
        p = ra.ranova.pValueGG_sel(ind); 
    end
    eta = ra.ranova.Eta2{ind};
    txt = sprintf('%s \n [F (%d,%d) = %.2f, p = %*.3f, %c2 = %0.2f]',ra.ranova.Row{ind},DF1,DF2,F,4,p,951,eta);
    ind = strfind(txt,'p = 0');
    txt(ind+4) = []; 
    ind = strfind(txt,'2 = 0');
    txt(ind+4) = []; 
    disp(txt)
    return;
end
try
    inds = find(ra.ranova.pValue_sel < 0.05);
catch
    inds = find(ra.ranova.pValueGG_sel < 0.05);
end
for ii = 1:length(inds)
    if inds(ii) == 1
        continue;
    end
    print_for_manuscript(ra,inds(ii));
end


function ind = row_identifier(rows,row)

for ii = 1:length(rows)
    thisr = rows{ii};
    indcol = strfind(rows{ii},':');
    if isempty(indcol)
        continue;
    end
    thisr = thisr((indcol(1)+1):end);
    if strcmp(thisr,row)
        break;
    end
end
ind = ii;
