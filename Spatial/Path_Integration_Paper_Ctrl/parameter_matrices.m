function [out,out1] = parameter_matrices(to_do,protocol,varargin)

p = inputParser;
addRequired(p,'to_do',@ischar);
addRequired(p,'protocol',@ischar);
addOptional(p,'data',[]);
parse(p,to_do,protocol,varargin{:});

data = p.Results.data;


out = [];
out1 = [];

if strcmp(lower(to_do),'calculate')
    get_parameters_matrices(protocol,data);
    return;
end

if strcmp(lower(to_do),'get')
    out = get_parameters_matrices(protocol);
    return;
end

if strcmp(lower(to_do),'select')
    [out,out1] = get_parameters_matrices(protocol,data{1},data{2});
    return;
end

if ~isempty(strfind(lower(to_do),'print'))
    cpMs = data{1};
    pMs = data{2}; 
    T = data{3};
    selAnimals = data{4};
    for rr = 1:size(pMs,1)
        for cc = 1:size(pMs,2)
            for ani = 1:length(selAnimals)
                an = selAnimals(ani);
                Perc_an(ani,rr,cc) = pMs{rr,cc}.perc(an);
                Num_an(ani,rr,cc) = pMs{rr,cc}.numCells(an);
            end
        end
    end
    perc_cells = squeeze(Perc_an);
    perc_cellsT = array2table(perc_cells);
    num_cells = squeeze(Num_an);
    num_cellsT = array2table(num_cells);
    if ~iscell(T{1,1})
        Ttext = table(num2str(T{selAnimals,1}),'VariableNames',{'Animal_id'});
    else
        Ttext = table((T{selAnimals,1}),'VariableNames',{'Animal_id'});
    end
    if strcmp(to_do,'print numbers')
        TT = [Ttext T(selAnimals,2) num_cellsT]
        out = num_cells;
    else
        TTp = [Ttext T(selAnimals,2) perc_cellsT]
        out = perc_cells;
    end
    ds = descriptiveStatistics(cpMs.perc(selAnimals));
    disp('Common');
    disp(cpMs.perc(selAnimals)');
    disp(cpMs.numCells(selAnimals)');
    disp('Average');
    disp(sprintf('Mean = %.3f - SEM = %.3f',ds.avg,ds.sem));
    return;
end