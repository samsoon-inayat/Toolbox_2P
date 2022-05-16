function txt = print_for_manuscript(ra,row)

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
    vartype = ra.ranova.Row{ind};
    ind = strfind(vartype,':');
    vtxt = vartype((ind(1)+1):end);
    inds = strfind(vtxt,'_');
    vtxt(inds) = '-';
    txt = sprintf('%s   [F (%d,%d) = %.2f, p = %*.3f, %c2 = %0.2f]',vtxt,DF1,DF2,F,4,p,951,eta);
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
if length(inds) > 1
    ht = 0.15*5;
    hf = get_figure(555,[10,10,3,ht]);
end
ys = 0; gap = 0.25;
for ii = 1:length(inds)
    if inds(ii) == 1
        continue;
    end
    txt = print_for_manuscript(ra,inds(ii));
    hdt = text(-0.1,ys,txt,'FontSize',6);
    te = get(hdt,'Extent');
    ys = te(2) + te(4) + gap;
end
axis off;
mData = evalin('base','mData');
save_pdf(hf,mData.pdf_folder,'rmanova.pdf',600);
% close(hf);
txt =[];

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
