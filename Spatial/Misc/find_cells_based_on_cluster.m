function [OK,classes,classesa,agfrmat] = find_cells_based_on_cluster(tresp)

good_FRmatD = double(tresp);

klist=2:10;%the number of clusters you want to try
myfunc = @(X,K)(kmeans(X, K));
eva = evalclusters(good_FRmatD,myfunc,'CalinskiHarabasz','klist',klist);
classes=kmeans(good_FRmatD,eva.OptimalK);

agfrmat = tresp(classes==1,:);
agfrmat = [agfrmat;tresp(classes==2,:)];

figure(200);clf;imagesc(agfrmat);colorbar;set(gca,'Ydir','normal')
OK = eva.OptimalK;