function temp_try
rng('default') % For reproducibility
X = rand(10,3);

tree = linkage(X,'average');
figure(10000);clf;
H = dendrogram(tree,'Orientation','left','ColorThreshold','default');
set(H,'LineWidth',2)