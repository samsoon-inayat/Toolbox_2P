function [pop, labels18] = build_pop18_singularB(idv, case_labels)
% Build 18 singular-only populations from idv (A×6 cell, each [nCells×3] logical [T D S]).
% A cell is in population (case t, label b) iff X(:,b)==1 AND sum(X,2)==1.
% Outputs:
%   pop.X{a,m} : logical vector (#cells(a)×1), m=1..18 in order [T D S] per case 1..6
%   labels18   : 1×18 labels like 'C3-A1-T', ... (or 'Case1-T' if not provided)

A = size(idv,1); C = size(idv,2); assert(C==4,'Expect 4 cases');
labs3 = {'T','M'};
if nargin<2 || isempty(case_labels)
    case_labels = arrayfun(@(k)sprintf('Case%d',k),1:C,'uni',0);
end

nCells = zeros(A,1);
for a=1:A, nCells(a) = size(idv{a,1},1); end
M = C*2; X = cell(A,M); labels18 = cell(1,M);

m = 0;
for t = 1:C
    for b = 1:2
        m = m + 1;
        labels18{m} = sprintf('%s-%s', case_labels{t}, labs3{b});
        for a = 1:A
            Xab = logical(idv{a,t});                   % [nCells×3]
            X{a,m} = (Xab(:,b)==1) & (sum(Xab,2)==1);  % singular-only
        end
    end
end

pop.A = A; pop.M = M; pop.nCells = nCells; pop.X = X;
end
