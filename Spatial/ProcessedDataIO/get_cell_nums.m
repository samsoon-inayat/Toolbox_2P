function cnso = get_cell_nums(Rs,resp,pl)

if isfield(Rs,'cns')
cns = Rs.cns;
c = cns(resp,:);
cnso = c(c(:,2)==pl,3);
else
    cnso = find(resp);
end