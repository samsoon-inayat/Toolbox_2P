function brake_vs_no_brake_cell_clustering


%%
good_FRmatD = double(good_FRmat);

[clus_inds,clus_centers] = kmeans(good_FRmat,4);
%%
klist=2:10;%the number of clusters you want to try
myfunc = @(X,K)(kmeans(X, K));
eva = evalclusters(good_FRmatD,myfunc,'CalinskiHarabasz','klist',klist)
classes=kmeans(good_FRmatD,eva.OptimalK);

agfrmat = good_FRmat(classes==1,:);
agfrmat = [agfrmat;good_FRmat(classes==2,:)];

figure(200);clf;imagesc(agfrmat);colorbar;set(gca,'Ydir','normal')