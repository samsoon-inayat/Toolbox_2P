function cluster_conj_comp(hf,mOI,txl,options)
mData = evalin('base','mData');
 %%
    mOI1 = mOI;
%     mOI1(isnan(mOI1)) = 1;
    Di = pdist(mOI1,@naneucdist);
%     tree = linkage(mOI1,'average','euclidean');
    tree = linkage(Di,'average');
    figure(hf);clf
    [H,T,TC] = dendrogram(tree,81,'Orientation','top','ColorThreshold','default');
    hf = gcf;
%     txl = event_type;
    set(hf,'Position',[7 3 6.9 1.5]);
    set(H,'linewidth',1);
    set(gca,'xticklabels',txl(TC));xtickangle(45);
    format_axes(gca);
    hx = ylabel('Eucledian Distance');changePosition(hx,[0 -0.1 0]);
    changePosition(gca,[-0.05 0.0 0.09 0.05]);