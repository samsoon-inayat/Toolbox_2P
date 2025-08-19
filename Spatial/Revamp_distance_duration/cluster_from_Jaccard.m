function [Z,ord,D] = cluster_from_Jaccard(Jmean)
% Jmean: 18×18 symmetric with diag≈1
D = 1 - Jmean; D(1:size(D,1)+1:end) = 0;   % zero diagonal
dvec = squareform(D);                      % condensed
Z = linkage(dvec,'average');
ord = optimalleaforder(Z, dvec);
end
