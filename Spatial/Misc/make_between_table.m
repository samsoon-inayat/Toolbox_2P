function between = make_between_table(data,dvn)
ngroups = size(data,1);

if ngroups == 1
    dataT = [];
    for cc = 1:size(data,2)
        dataT = [dataT data{cc}];
    end
    between = array2table(dataT);
    between.Properties.VariableNames = dvn;
    return;
end

if ngroups > 1
    dataT = [];
    group = [];
    for rr = 1:size(data,1)
        dataC = [];
        for cc = 1:size(data,2)
            dataC = [dataC data{rr,cc}];
        end
        dataT = [dataT;dataC];
        group = [group;(ones(size(dataC,1),1)*rr)];
    end
    between = array2table([group dataT]);
    between.Properties.VariableNames = {'Group',dvn{:}};
    between.Group = categorical(between.Group);
end

