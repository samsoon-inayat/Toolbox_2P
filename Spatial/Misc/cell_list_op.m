function resp = cell_list_op(resp_valsCi,resp_valsC1i,fun)

if isstruct(resp_valsCi)
    resp_valsC = resp_valsCi.vals;
    resp_valsC1 = resp_valsC1i.vals;
else
    resp_valsC = resp_valsCi;
    resp_valsC1 = resp_valsC1i;
end

if ~isempty(strfind(fun,'not'))
    for rr = 1:size(resp_valsC,1)
        for cc = 1:size(resp_valsC,2)
            resp{rr,cc} = ~resp_valsC{rr,cc};
        end
    end
    return;
end

if ~isempty(strfind(fun,'only'))
    ind = str2num(fun(5:end));
    all = 1:size(resp_valsC,2);
    sd = setdiff(all,ind);
    resp_sd = cell_list_op(cell_list_op(resp_valsC(:,sd),[],'or'),[],'not');
    resp = cell_list_op(resp_valsC(:,1),resp_sd(:,1),'and');
    resp = repmat(resp,1,size(resp_valsC,2));
    return;
end


if ~isempty(resp_valsC1i)
    for rr = 1:size(resp_valsC,1)
        for cc = 1:size(resp_valsC,2)
            c = resp_valsC{rr,cc}; c1 = resp_valsC1{rr,cc};
            if strcmp(fun,'and')
                resp{rr,cc} = c & c1;
            end
            if strcmp(fun,'or')
                resp{rr,cc} = c | c1;
            end
        end
    end
else
    for rr = 1:size(resp_valsC,1)
        if strcmp(fun,'or')
            respTemp = logical(zeros(size(resp_valsC{rr,1})));
        end
        if strcmp(fun,'and')
            respTemp = logical(ones(size(resp_valsC{rr,1})));
        end
        for cc = 1:size(resp_valsC,2)
            if strcmp(fun,'or')
                respTemp = respTemp | resp_valsC{rr,cc};
            end
            if strcmp(fun,'and')
                respTemp = respTemp & resp_valsC{rr,cc};
            end
        end
        resp{rr} = respTemp;
    end
    resp = repmat(resp',1,size(resp_valsC,2));
end