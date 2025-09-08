function S = pop18_agreement_cluster(pop)
% 18×18 Agreement (= Jaccard) per animal; group mean/SEM; clustering on 1−Agreement.
A = pop.A; M = pop.M;
J_animal = nan(M,M,A); N_union = nan(M,M,A);

for a = 1:A
    n = pop.nCells(a);
    % stack membership to [n × M]
    Xmat = false(n,M); for m=1:M, xi = pop.X{a,m}; Xmat(1:numel(xi),m) = xi(:); end
    J = nan(M,M); Uc = nan(M,M);
    for i=1:M
        Ai = Xmat(:,i);
        for j=1:M
            if i==j, J(i,j)=NaN; Uc(i,j)=NaN; continue; end
            Aj = Xmat(:,j);
            k = sum(Ai & Aj);
            u = sum(Ai) + sum(Aj) - k;
            if u==0, J(i,j)=NaN; Uc(i,j)=0; else, J(i,j)=k/u; Uc(i,j)=u; end
        end
    end
    J_animal(:,:,a) = J; N_union(:,:,a) = Uc;
end

J_mean = mean(J_animal,3,'omitnan');
na     = sum(~isnan(J_animal),3);
J_sem  = std(J_animal,0,3,'omitnan') ./ sqrt(max(na,1));

% clustering on Distance = 1 − Agreement
D = 1 - J_mean; D(1:M+1:end)=0; D(isnan(D))=0;
Y = squareform(D,'tovector');
tree = linkage(Y,'average');
leafOrder = optimalleaforder(tree,Y);
[coph_r,~] = cophenet(tree,Y);

S.J_animal = J_animal; S.J_mean = J_mean; S.J_sem = J_sem;
S.N_union = N_union; S.tree = tree; S.leafOrder = leafOrder; S.coph_r = coph_r;
end
