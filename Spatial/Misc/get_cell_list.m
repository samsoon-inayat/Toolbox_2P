function [selResp,selResp2] = get_cell_list(resp_valsCi,cond,per)

% resp_valsCi = resp_valsC;
for rr = 1:size(resp_valsCi,1)
    ccs = [];
    for cc = 1:size(resp_valsCi,2)
        ccs(:,cc) = resp_valsCi{rr,cc};
    end
    resp_valsC{rr} = ccs;
end


for rr = 1:size(resp_valsCi,1)
    for cc = 1:size(resp_valsCi,2)
        selResp2{rr,cc} = resp_valsCi{rr,cc};
    end
end


if isempty(cond)
    for ii = 1:length(resp_valsC)
        selResp{ii} = logical(ones(size(resp_valsC{ii},1),1));
    end
    return;
end

for ii = 1:length(resp_valsC)
    selResp{ii} = logical(zeros(size(resp_valsC{ii},1),1)); % final response result holder
    for rr = 1:size(cond,1)
        temp_resp = logical(ones(size(resp_valsC{ii},1),1)); % make all logical ones
        for cc = 1:size(cond,2) % if there are cols in the cond then they will be ANDed together.
            if isnan(cond(rr,cc)) % if cond is NaN disregard that condition
                this_list = logical(ones(size(resp_valsC{ii},1),1));
            end
            if cond(rr,cc) < 0 % if cond is negative means do the not operation
                this_list = ~logical(resp_valsC{ii}(:,abs(cond(rr,cc))));
            end
            if cond(rr,cc) > 0
                this_list = logical(resp_valsC{ii}(:,abs(cond(rr,cc))));
            end
            temp_resp = temp_resp & this_list;
        end
        selResp{ii} = selResp{ii} | temp_resp; % multiple rows in cond will be ORed together
    end
end

if exist('per','var')
    if per == 0
        selResp = selResp';
        return;
    end
    for rr = 1:length(selResp)
        perc(rr) = 100*sum(selResp{rr})/length(selResp{rr});
    end
    selResp = perc;
end

selResp = selResp';
selResp = repmat(selResp,1,size(resp_valsCi,2));
