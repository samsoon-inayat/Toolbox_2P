function selResp = get_cell_list(resp_valsC,cond)

if isempty(cond)
    for ii = 1:length(resp_valsC)
        selResp{ii} = logical(ones(size(resp_valsC{ii},1),1));
    end
    return;
end

for ii = 1:length(resp_valsC)
    selResp{ii} = logical(zeros(size(resp_valsC{ii},1),1));
    for rr = 1:size(cond,1)
        temp_resp = logical(ones(size(resp_valsC{ii},1),1));
        for cc = 1:size(cond,2)
            if isnan(cond(rr,cc))
                this_list = logical(ones(size(resp_valsC{ii},1),1));
            end
            if cond(rr,cc) < 0
                this_list = ~logical(resp_valsC{ii}(:,abs(cond(rr,cc))));
            end
            if cond(rr,cc) > 0
                this_list = logical(resp_valsC{ii}(:,abs(cond(rr,cc))));
            end
            temp_resp = temp_resp & this_list;
        end
        selResp{ii} = selResp{ii} | temp_resp;
    end
end