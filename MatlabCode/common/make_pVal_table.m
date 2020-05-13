function hf = make_pVal_table(c,legendT_all,figN,figPos)

for ii = 1:size(c,1)
    this_p_val = c(ii,end);
    p_str = getNumberOfAsterisks(this_p_val);
    pValTab{c(ii,1),c(ii,2)} = sprintf('%s',p_str);
%     pValTab{c(ii,1),c(ii,2)} = sprintf('%s (%.3f)',p_str,this_p_val);
end

hf = figure(figN);
% set(hf,'Position',[50 50 800 500]);
set(hf,'Position',figPos);
tbl = uitable(hf,'Position',[20 20 800 300]);
tbl.Data = pValTab;
tbl.RowName = legendT_all;tbl.ColumnName = legendT_all;
% for ii = 1:7
%     tbl.RowName{ii} = sprintf('<html><font size=+1><b>%s</b></font></html>',legendT_all{ii});
%     tbl.ColumnName{ii} = sprintf('<html><font size=+1><b>%s</b></font></html>',legendT_all{ii});
% end
tbl.FontSize = 10;