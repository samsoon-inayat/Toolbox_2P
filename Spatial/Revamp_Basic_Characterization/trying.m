function trying

%% categorize cells into clusters from trial responses
while 1
    
    tresp = cell2mat(allresp(1,1:10));
    figure(100);clf;imagesc(tresp)
    good_FRmatD = double(tresp);

    [clus_inds,clus_centers] = kmeans(tresp,4);

    klist=2:10;%the number of clusters you want to try
    myfunc = @(X,K)(kmeans(X, K));
    eva = evalclusters(good_FRmatD,myfunc,'CalinskiHarabasz','klist',klist)
    classes=kmeans(good_FRmatD,eva.OptimalK);

    agfrmat = tresp(classes==1,:);
    agfrmat = [agfrmat;tresp(classes==2,:)];

    figure(200);clf;imagesc(agfrmat);colorbar;set(gca,'Ydir','normal')
    100*sum(classes==2)/size(tresp,1)
    100*sum(classes==1)/size(tresp,1)
    break;
end

%% testing whether the method of clustering is fine for finding responsive cells
sz1 = size(trialR,1);
sz2 = size(trialR,2);
trialRa = [trialR(FR_based == 0,:);trialR(FR_based == 1,:)];
fd = BoxCountfracDim(trialRa);
for ii = 1:10
    trialRS = [];
    for jj = 1:sz2
        trialRS = [trialRS trialRa(randperm(sz1,sz1),jj)];
    end
    ocs(ii) = find_cells_based_on_cluster(trialRS);
    fds(ii) = BoxCountfracDim(trialRS);
end
[h,p,ci,stats] = ttest(fds,fd);
[h p]
%%
tR = rasters(:,:,7)';
trialR = tR>0;
oc = 2;
rng(1);
sRp = kmeans(trialR,2);
s1 = median(sum(trialR(sRp==1,:),2));
s2 = median(sum(trialR(sRp==2,:),2));
if s1>s2
    sRp1 = sRp == 1;
else
    sRp1 = sRp == 2;
end
FR_based = sRp1;

%%
for fi = 1:500%length(frames)
    figure(100);clf;
    imagesc(frames{fi});
    title(fi);
    pause(0.1);
end

