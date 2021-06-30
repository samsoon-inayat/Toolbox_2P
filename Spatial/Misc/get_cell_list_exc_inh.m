function [resp,resp_OR,resp_AND,resp_fraction] = get_cell_list_exc_inh(resp_exc_inh,cond,exc_inh)

for ii = 1:length(resp_exc_inh)
    tresp1 = resp_exc_inh{ii,1}(:,1);
    resp_OR{ii} = logical(zeros(size(tresp1,1),1));
    resp_AND{ii} = logical(ones(size(tresp1,1),1));
    for cc = 1:size(resp_exc_inh,2)
        if exc_inh
            tresp = resp_exc_inh{ii,1}(:,cc);
        else
            tresp = resp_exc_inh{ii,2}(:,cc);
        end
        resp_OR{ii} = resp_OR{ii} | tresp;
        resp_AND{ii} = resp_AND{ii} & tresp;
        if cc == cond
            resp{ii} = tresp;
            resp_fraction(ii) = sum(tresp)/length(tresp);
        end
    end
    n = 0;
end