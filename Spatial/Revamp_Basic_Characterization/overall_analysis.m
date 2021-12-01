function overall_analysis

%%
si = [Lb_T ArL_L_T Lbs_T Ab_T Abs_T Ar_t_D ArL_t_D Ars_t_D Ar_i_T ArL_i_T Ars_i_T];
% si = [Ab_T Abs_T];
% si = [Ar_t_D ArL_t_D Ars_t_D ];
Rs = o.Rs(:,si); mRs = o.mR(:,si);
% 
% si = [Ab_T Abs_T Ar_t_D ArL_t_D Ars_t_D Ar_i_T ArL_i_T Ars_i_T];
% RsR = o.Rs(:,si);
props1 = get_props_Rs(Rs,50);

% resp_all = props1.good_FR;
% resp = get_cell_list(resp_all,[3]);
% resp = [resp(:,1:3) resp];
resp = props1.good_Gauss_loose;
view_population_vector(Rs,mRs,resp,100);
%%
view_population_vector_corr(Rs,mRs,1,200);

%% Overlap Indices ImageSC
while 1
    ntrials = 50;
    fileName = fullfile(mData.pd_folder,sprintf('%s_tuned_cells',mfilename));
    load(fileName);
    si = [MOn_T MOff_T Lb_T ArL_L_T Lbs_T Ab_T Abs_T Ar_t_D ArL_t_D Ars_t_D Ar_i_T ArL_i_T Ars_i_T];
    si = [MOn_T MOff_T Lb_T ArL_L_T Lbs_T Ab_T Abs_T Ar_t_D ArL_t_D Ars_t_D Ar_i_T ArL_i_T Ars_i_T];
    props1 = get_props_Rs(o.Rs(:,si),ntrials);
    resp = [ speed_tuned_cells props1.good_FR];% resp_speed];
    si = [Ab_t_T Ab_i_T Abs_t_T Abs_i_T]; 
    props1 = get_props_Rs(o.Rs(:,si),ntrials);
    resp = props1.good_FR;
   
%     resp(:,6:8) = dzMI.resp_D_g_T(:,[1 3 5]); resp(:,9:11) = dzMI.resp_T_g_D(:,[2 4 6]);
%     resp(:,6:8) = dzMI.resp_D_g_T(:,[1 3 5]); resp(:,9:11) = dzMI.resp_D_g_T(:,[2 4 6]);
%     resp(:,6:8) = dzMI.resp_T_g_D(:,[1 3 5]); resp(:,9:11) = dzMI.resp_T_g_D(:,[2 4 6]);
    [OI,mOI,semOI,OI_mat,p_vals,h_vals] = get_overlap_index(resp,0.5,0.05);
    sz = size(mOI,1);
%     mOI = OI_mat(:,:,4);
    oM = ones(size(mOI));
    mask1 = (triu(oM,0) & tril(oM,0)); mOI(mask1==1) = NaN;
    maxI = max([mOI(:);semOI(:)]);    
    minI = min([mOI(:);semOI(:)]);
    
    mask = tril(NaN(size(mOI)),0); mask(mask==0) = 1; 
    txl = [{'sR'} rasterNamesTxt(si)]; 
    txl = [rasterNamesTxt(si)]; 
%     mOI = mOI .* mask;
    imAlpha=ones(size(mOI));    %imAlpha(isnan(mask))=0.25; 
    imAlpha(mask1 == 1) = 0;
%     ff = makeFigureRowsCols(2020,[10 4 6 1.5],'RowsCols',[1 2],'spaceRowsCols',[0.1 0.01],'rightUpShifts',[0.05 0.13],'widthHeightAdjustment',[-240 -150]);
    hf = get_figure(5,[8 7 1.5 1.5]);
    %
%     axes(ff.h_axes(1));
    im1 = imagesc(mOI,[minI,maxI]);    im1.AlphaData = imAlpha;
    set(gca,'color',0.5*[1 1 1]);    colormap parula;    %axis equal
    format_axes(gca);
    set_axes_limits(gca,[0.5 sz+0.5],[0.5 sz+0.5]);
    set(gca,'xtick',1:length(txl),'ytick',1:length(txl),'xticklabels',txl,'yticklabels',txl,'Ydir','reverse'); xtickangle(45);
    text(0.5,sz+1.1,'Average Overlap Index (N = 5 animals)','FontSize',5); set(gca,'Ydir','normal');ytickangle(20);
    box on
    changePosition(gca,[0.03 0 -0.04 0]);
    hc = putColorBar(gca,[0.09 0.07 -0.11 -0.15],{sprintf('%d',minI),sprintf('%.1f',maxI)},6,'eastoutside',[0.1 0.05 0.06 0.05]);
    save_pdf(hf,mData.pdf_folder,sprintf('OI_Map_%d_mean.pdf',ntrials),600);
    %%
    break;
end

%% agglomerative hierarchical clustering
while 1
    mOI1 = mOI;
    mOI1(isnan(mOI1)) = 1;
    Di = pdist(mOI1);
    tree = linkage(Di);
    figure(hf);clf
    [H,T,TC] = dendrogram(tree,'Orientation','right','ColorThreshold','default');
    hf = gcf;
    set(hf,'Position',[7 3 1.25 2]);
    set(H,'linewidth',1);
    set(gca,'yticklabels',txl(TC));ytickangle(30);
    format_axes(gca);
    hx = xlabel('Eucledian Distance');%changePosition(hx,[-0.051 0 0]);
    changePosition(gca,[0 0.0 0.05 0.05]);
    save_pdf(hf,mData.pdf_folder,sprintf('OI_Map_cluster.pdf'),600);
    %%
    break;
end



%% Overlap Indices ImageSC check for overlap among air stream rest versus motion trials and intertrials
while 1
    ntrials = 50;
%     fileName = fullfile(mData.pd_folder,sprintf('%s_tuned_cells',mfilename));
%     load(fileName);
    
    si = [Ab_t_T Abs_t_T Ar_t_D ArL_t_D Ars_t_D Ab_i_T Abs_i_T Ar_i_T ArL_i_T Ars_i_T];
    props1 = get_props_Rs(o.Rs(:,si),ntrials);
    
    resp = [props1.good_FR];% resp_speed];
%     resp(:,6:8) = dzMI.resp_D_g_T(:,[1 3 5]); resp(:,9:11) = dzMI.resp_T_g_D(:,[2 4 6]);
%     resp(:,6:8) = dzMI.resp_D_g_T(:,[1 3 5]); resp(:,9:11) = dzMI.resp_D_g_T(:,[2 4 6]);
%     resp(:,6:8) = dzMI.resp_T_g_D(:,[1 3 5]); resp(:,9:11) = dzMI.resp_T_g_D(:,[2 4 6]);
    [OI,mOI,semOI,OI_mat,p_vals,h_vals] = get_overlap_index(resp,0.5,0.05);
    sz = size(mOI,1);
%     mOI = OI_mat(:,:,4);
    oM = ones(size(mOI));
    mask1 = (triu(oM,0) & tril(oM,0)); mOI(mask1==1) = NaN;
    maxI = max([mOI(:);semOI(:)]);    
    minI = min([mOI(:);semOI(:)]);
    
    mask = tril(NaN(size(mOI)),0); mask(mask==0) = 1; 
    txl = [rasterNamesTxt(si)]; 
%     mOI = mOI .* mask;
    imAlpha=ones(size(mOI));    %imAlpha(isnan(mask))=0.25; 
    imAlpha(mask1 == 1) = 0;
%     ff = makeFigureRowsCols(2020,[10 4 6 1.5],'RowsCols',[1 2],'spaceRowsCols',[0.1 0.01],'rightUpShifts',[0.05 0.13],'widthHeightAdjustment',[-240 -150]);
    hf = get_figure(5,[8 7 1.5 1.5]);
    %
%     axes(ff.h_axes(1));
    im1 = imagesc(mOI,[minI,maxI]);    im1.AlphaData = imAlpha;
    set(gca,'color',0.5*[1 1 1]);    colormap parula;    %axis equal
    format_axes(gca);
    set_axes_limits(gca,[0.5 sz+0.5],[0.5 sz+0.5]);
    set(gca,'xtick',1:length(txl),'ytick',1:length(txl),'xticklabels',txl,'yticklabels',txl,'Ydir','reverse'); xtickangle(45);
    text(0.5,sz+1.1,'Average Overlap Index (N = 5 animals)','FontSize',5); set(gca,'Ydir','normal');ytickangle(20);
    box on
    changePosition(gca,[0.03 0 -0.04 0]);
    hc = putColorBar(gca,[0.09 0.07 -0.11 -0.15],{sprintf('%d',minI),sprintf('%.1f',maxI)},6,'eastoutside',[0.1 0.05 0.06 0.05]);
    save_pdf(hf,mData.pdf_folder,sprintf('OI_Map_%d_mean.pdf',ntrials),600);
    %%
    break;
end

%% agglomerative hierarchical clustering
while 1
    mOI1 = mOI;
    mOI1(isnan(mOI1)) = 1;
    Di = pdist(mOI1);
    tree = linkage(Di);
    figure(hf);clf
    [H,T,TC] = dendrogram(tree,'Orientation','right','ColorThreshold','default');
    hf = gcf;
    set(hf,'Position',[7 3 1.25 2]);
    set(H,'linewidth',1);
    set(gca,'yticklabels',txl(TC));ytickangle(30);
    format_axes(gca);
    hx = xlabel('Eucledian Distance');%changePosition(hx,[-0.051 0 0]);
    changePosition(gca,[0 0.0 0.05 0.05]);
    save_pdf(hf,mData.pdf_folder,sprintf('OI_Map_cluster.pdf'),600);
    %%
    break;
end


%% Overlap Indices ImageSC
while 1
    ntrials = 50;
    fileName = fullfile(mData.pd_folder,sprintf('%s_tuned_cells',mfilename));
    load(fileName);
    si = [MOn_T MOff_T Lb_T ArL_L_T Lbs_T Ab_T Abs_T Ar_t_D ArL_t_D Ars_t_D Ar_i_T ArL_i_T Ars_i_T];
    props1 = get_props_Rs(o.Rs(:,si),ntrials);
    resp = [ speed_tuned_cells props1.good_FR];% resp_speed];
   
%     resp(:,6:8) = dzMI.resp_D_g_T(:,[1 3 5]); resp(:,9:11) = dzMI.resp_T_g_D(:,[2 4 6]);
%     resp(:,6:8) = dzMI.resp_D_g_T(:,[1 3 5]); resp(:,9:11) = dzMI.resp_D_g_T(:,[2 4 6]);
%     resp(:,6:8) = dzMI.resp_T_g_D(:,[1 3 5]); resp(:,9:11) = dzMI.resp_T_g_D(:,[2 4 6]);
    [OI,mOI,semOI,OI_mat,p_vals,h_vals] = get_overlap_index(resp,0.5,0.05);
    sz = size(mOI,1);
%     mOI = OI_mat(:,:,4);
    oM = ones(size(mOI));
    mask1 = (triu(oM,0) & tril(oM,0)); mOI(mask1==1) = NaN;
    maxI = max([mOI(:);semOI(:)]);    
    minI = min([mOI(:);semOI(:)]);
    
    mask = tril(NaN(size(mOI)),0); mask(mask==0) = 1; 
    txl = [{'sR'} rasterNamesTxt(si)]; 
%     mOI = mOI .* mask;
    imAlpha=ones(size(mOI));    %imAlpha(isnan(mask))=0.25; 
    imAlpha(mask1 == 1) = 0;
%     ff = makeFigureRowsCols(2020,[10 4 6 1.5],'RowsCols',[1 2],'spaceRowsCols',[0.1 0.01],'rightUpShifts',[0.05 0.13],'widthHeightAdjustment',[-240 -150]);
    hf = get_figure(5,[8 7 1.5 1.5]);
    %
%     axes(ff.h_axes(1));
    im1 = imagesc(mOI,[minI,maxI]);    im1.AlphaData = imAlpha;
    set(gca,'color',0.5*[1 1 1]);    colormap parula;    %axis equal
    format_axes(gca);
    set_axes_limits(gca,[0.5 sz+0.5],[0.5 sz+0.5]);
    set(gca,'xtick',1:length(txl),'ytick',1:length(txl),'xticklabels',txl,'yticklabels',txl,'Ydir','reverse'); xtickangle(45);
    text(0.5,sz+1.1,'Average Overlap Index (N = 5 animals)','FontSize',5); set(gca,'Ydir','normal');ytickangle(20);
    box on
    changePosition(gca,[0.03 0 -0.04 0]);
    hc = putColorBar(gca,[0.09 0.07 -0.11 -0.15],{sprintf('%d',minI),sprintf('%.1f',maxI)},6,'eastoutside',[0.1 0.05 0.06 0.05]);
    save_pdf(hf,mData.pdf_folder,sprintf('OI_Map_%d_mean.pdf',ntrials),600);
    %%
    break;
end

%% agglomerative hierarchical clustering
while 1
    mOI1 = mOI;
    mOI1(isnan(mOI1)) = 1;
    Di = pdist(mOI1);
    tree = linkage(Di);
    figure(hf);clf
    [H,T,TC] = dendrogram(tree,'Orientation','right','ColorThreshold','default');
    hf = gcf;
    set(hf,'Position',[7 3 1.25 2]);
    set(H,'linewidth',1);
    set(gca,'yticklabels',txl(TC));ytickangle(30);
    format_axes(gca);
    hx = xlabel('Eucledian Distance');%changePosition(hx,[-0.051 0 0]);
    changePosition(gca,[0 0.0 0.05 0.05]);
    save_pdf(hf,mData.pdf_folder,sprintf('OI_Map_cluster.pdf'),600);
    %%
    break;
end
