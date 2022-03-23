function distance_duration_pops

%%
ntrials = [50];
RsDt = o.Rs(:,[Ar_t_D ArL_t_D Ars_t_D]);  RsTt = o.Rs(:,[Ar_t_T ArL_t_T Ars_t_T]);
RsDi = o.Rs(:,[Ar_i_D ArL_i_D Ars_i_D]);  RsTi = o.Rs(:,[Ar_i_T ArL_i_T Ars_i_T]);
[dzMI_FD,dzMI_FT] = get_zMI_comp_dist_time(RsDt,RsTt,RsDi,RsTi);

[allRs_FD,allmRs_FD] = get_trial_Rs(o,[Ar_t_D ArL_t_D Ars_t_D],1:10);
[allRs_FT,allmRs_FT] = get_trial_Rs(o,[Ar_i_T ArL_i_T Ars_i_T],1:10);
[allRs_Ab,allmRs_Ab] = get_trial_Rs(o,[Ab_T Abs_T],1:10);
respDTC = combine_distance_time_rasters(o.Rs(:,[Ar_t_D ArL_t_D Ars_t_D]),o.Rs(:,[Ar_i_T ArL_i_T Ars_i_T]),ntrials);


mRsDt = o.mR(:,[Ar_t_D ArL_t_D Ars_t_D]);  mRsTt = o.mR(:,[Ar_t_T ArL_t_T Ars_t_T]);
mRsDi = o.mR(:,[Ar_i_D ArL_i_D Ars_i_D]);  mRsTi = o.mR(:,[Ar_i_T ArL_i_T Ars_i_T]);

Rs_Ab = o.Rs(:,[Ab_T Abs_T]); mRs_Ab = o.mR(:,[Ab_T Abs_T]);
%%
trialsR = [50];
props_Ab = get_props_Rs(Rs_Ab,trialsR); propsD = get_props_Rs(RsDt,trialsR); propsT = get_props_Rs(RsTi,trialsR);

trialsR = [10 30];
props_Ab13 = get_props_Rs(Rs_Ab,trialsR); propsD13 = get_props_Rs(RsDt,trialsR); propsT13 = get_props_Rs(RsTi,trialsR);

trialsR = [40 70];
props_Ab47 = get_props_Rs(Rs_Ab,trialsR); propsD47 = get_props_Rs(RsDt,trialsR); propsT47 = get_props_Rs(RsTi,trialsR);

trialsR = [80 100];
props_Ab810 = get_props_Rs(Rs_Ab,trialsR); propsD810 = get_props_Rs(RsDt,trialsR); propsT810 = get_props_Rs(RsTi,trialsR);

%% visualizing different cell groups for the same condition
while 1
    %%
    an = 3; cn = 1;
    
    tRs = [RsDt(an,cn) allRs_FD{cn}(an,:)];  tmRs = [mRsDt(an,cn) allmRs_FD{cn}(an,:)];
    tRsi = [RsTi(an,cn) allRs_FT{cn}(an,:)];  tmRsi = [mRsTi(an,cn) allmRs_FT{cn}(an,:)];

    resp = respDT.inh(an,cn);
%     resp = respDT.exc_and_good_FR(an,cn);
    resp = propsD13.good_FR(an,cn);
    ff = plot_pop_vec_trials(100,[1 1 14 3],tRs,tmRs,resp)
%     colormap parula
    
    resp = respDT.inh_and_good_FR(an,cn);
%     resp = respDT.exc_not_good_FR(an,cn);
    resp = propsD47.good_FR(an,cn);
    ff = plot_pop_vec_trials(101,[1 4 14 3],tRsi,tmRsi,resp)

    
    resp = propsD810.good_FR(an,cn);
    ff = plot_pop_vec_trials(102,[1 7 14 3],tRsi,tmRsi,resp)

    %%
    break;
end

%%
while 1
    %%
    cns = 3;
    an = 1; cn = cns(1);
    cnab = 2;
    resp = respDT.exc(an,cn);
    resp = props_Ab.good_FR(an,cnab);
    resp = propsD.good_FR(an,cn);
%     resp = dzMI_FD.resp_D_g_T_and_good_FR(an,cn);
    
    
    tRs = [RsDt(an,cn) allRs_FD{cn}(an,:)];  tmRs = [mRsDt(an,cn) allmRs_FD{cn}(an,:)];
    ff = plot_pop_vec_trials(100,[1 1 14 3],tRs,tmRs,resp)
    colormap jet
    
    tRs = [RsTi(an,cn) allRs_FT{cn}(an,:)];  tmRs = [mRsTi(an,cn) allmRs_FT{cn}(an,:)];
    ff = plot_pop_vec_trials(101,[1 4 14 3],tRs,tmRs,resp)
    colormap jet
    
    tRs = [Rs_Ab(an,cnab) allRs_Ab{cnab}(an,:)];  tmRs = [mRs_Ab(an,cnab) allmRs_Ab{cnab}(an,:)];
%     resp = respDT.exc(an,cn);
    ff = plot_pop_vec_trials(102,[1 7 14 3],tRs,tmRs,resp)
    colormap jet
    %%
    break;
end;

%% Overlap Indices ImageSC overall
while 1
    cn = 1;
    si = [Ar_t_D ArL_t_D Ars_t_D Ar_i_T ArL_i_T Ars_i_T Ar_t_D ArL_t_D Ars_t_D Ar_i_T ArL_i_T Ars_i_T];
    resp  = [respDT.inh respDT.exc propsD.good_FR propsT.good_FR]
    resp  = [respDT.inh respDT.exc]
    resp = [props_Ab.good_FR(:,cn),propsD.good_FR(:,cn),propsT.good_FR(:,cn)];
%     resp = [props_Ab13.good_FR,propsD13.good_FR,propsT13.good_FR];
    [OIo,mOI,semOI,OI_mato,p_vals,h_vals,all_CI,mCI,semCI,all_CI_mat,uni] = get_overlap_index(resp,0.5,0.05);
    mOI = mCI;
    sz = size(mOI,1);
%     mOI = OI_mat(:,:,4);
    oM = ones(size(mOI));
    mask1 = (triu(oM,0) & tril(oM,0)); mOI(mask1==1) = NaN;
    maxI = max([mOI(:)]);%semOI(:)]);    
    minI = min([mOI(:)]);%semOI(:)]);
    
    mask = tril(NaN(size(mOI)),0); mask(mask==0) = 1;
%     ptxl = {'D','D','D','T','T','T','C','C','C','D','D','D','T','T','T','C','C','C'};
    ptxl = {'D','D','D','T','T','T','','','','','',''};
    ptxl = {'D','D','D','','',''};
%     ptxl = {'D','T','C','D','T','C'};
    for ii = 1:length(ptxl)
        txl{ii} = sprintf('%s-%s',ptxl{ii},rasterNamesTxt{si(ii)});
    end
    txl = rasterNamesTxt([Ab_T Abs_T Ar_t_D ArL_t_D Ars_t_D Ar_i_T ArL_i_T Ars_i_T]);
    txl = rasterNamesTxt([Ab_T Ar_t_D Ar_i_T]);
%     txl = {'1-13','2-13','3-13','1-47','2-47','3-47','1-81','2-81','3-81'};
   
%     mOI = mOI .* mask;
    imAlpha=ones(size(mOI));    %imAlpha(isnan(mask))=0.25; 
    imAlpha(mask1 == 1) = 0;
%     ff = makeFigureRowsCols(2020,[10 4 6 1.5],'RowsCols',[1 2],'spaceRowsCols',[0.1 0.01],'rightUpShifts',[0.05 0.13],'widthHeightAdjustment',[-240 -150]);
    hf = get_figure(5,[8 3 3.5 3.5]);
    %
%     axes(ff.h_axes(1));
    im1 = imagesc(mOI,[minI,maxI]);    im1.AlphaData = imAlpha;
    set(gca,'color',0.5*[1 1 1]);    colormap parula;    %axis equal
    format_axes(gca);
    for rr = 1:size(mCI,1)
        for cc = 1:size(mCI,1)
            if rr == cc
                text(rr-0.25,cc,sprintf('%.0f',mCI(rr,cc)),'FontSize',5);
            end
        end
    end
    axis equal
    set_axes_limits(gca,[0.5 sz+0.5],[0.5 sz+0.5]);
%     set(gca,'xtick',1:length(txl),'ytick',1:length(txl),'xticklabels',txl,'yticklabels',txl,'Ydir','reverse'); xtickangle(45);
    set(gca,'xtick',1:length(txl),'ytick',1:length(txl),'xticklabels',txl,'yticklabels',txl,'Ydir','reverse'); xtickangle(45);
%     text(-0.5,sz+1.1,'Average Overlap Index (N = 5 animals)','FontSize',5); 
    set(gca,'Ydir','normal');ytickangle(25);
    box on
    changePosition(gca,[0.0 0 -0.04 0]);
    hc = putColorBar(gca,[0.08 0.07 -0.11 -0.15],{sprintf('%.1f',minI),sprintf('%.1f',maxI)},6,'eastoutside',[0.1 0.08 0.11 0.15]);
    colormap jet
    save_pdf(hf,mData.pdf_folder,sprintf('OI_Map_%d_mean1_conjunctionP.pdf',ntrials),600);
    break;
end

%% Overlap Indices ImageSC
while 1
    sii = [Ar_t_D ArL_t_D Ars_t_D Ar_t_D ArL_t_D Ars_t_D Ar_t_D ArL_t_D Ars_t_D Ar_i_T ArL_i_T Ars_i_T Ar_i_T ArL_i_T Ars_i_T Ar_i_T ArL_i_T Ars_i_T];
    sii = [Ar_t_D Ar_t_D Ar_t_D Ar_i_T Ar_i_T Ar_i_T];
    coln = 1;
    resp = [FDgFR_D_g_T(:,coln) FDgFR_T_g_D(:,coln) FDgFR_Comp(:,coln) FTgFR_D_g_T(:,coln) FTgFR_T_g_D(:,coln) FTgFR_Comp(:,coln)];

%     [OI,mOI,semOI,OI_mat,p_vals,h_vals] = get_overlap_index(resp,0.5,0.05);
    [OIo,mOI,semOI,OI_mato,p_vals,h_vals,all_CI,mCI,semCI,all_CI_mat,uni] = get_overlap_index(resp,0.5,0.05);
    mOI = mCI;
    sz = size(mOI,1);
%     mOI = OI_mat(:,:,4);
    oM = ones(size(mOI));
    mask1 = (triu(oM,0) & tril(oM,0)); mOI(mask1==1) = NaN;
    maxI = max([mOI(:)]);%semOI(:)]);    
    minI = min([mOI(:)]);%semOI(:)]);
    
    mask = tril(NaN(size(mOI)),0); mask(mask==0) = 1;
%     ptxl = {'D','D','D','T','T','T','C','C','C','D','D','D','T','T','T','C','C','C'};
    ptxl = {'D','T','C','D','T','C'};
    for ii = 1:length(ptxl)
        txl{ii} = sprintf('%s-%s',ptxl{ii},rasterNamesTxt{sii(ii)});
    end

    
%     mOI = mOI .* mask;
    imAlpha=ones(size(mOI));    %imAlpha(isnan(mask))=0.25; 
    imAlpha(mask1 == 1) = 0;
%     ff = makeFigureRowsCols(2020,[10 4 6 1.5],'RowsCols',[1 2],'spaceRowsCols',[0.1 0.01],'rightUpShifts',[0.05 0.13],'widthHeightAdjustment',[-240 -150]);
    hf = get_figure(5,[8 3 3.5 3.5]);
    %
%     axes(ff.h_axes(1));
    im1 = imagesc(mOI,[minI,maxI]);    im1.AlphaData = imAlpha;
    set(gca,'color',0.5*[1 1 1]);    colormap parula;    %axis equal
    format_axes(gca);
    for rr = 1:size(mCI,1)
        for cc = 1:size(mCI,1)
            if rr == cc
                text(rr-0.25,cc,sprintf('%.0f',mCI(rr,cc)),'FontSize',5);
            end
        end
    end
    axis equal
    set_axes_limits(gca,[0.5 sz+0.5],[0.5 sz+0.5]);
%     set(gca,'xtick',1:length(txl),'ytick',1:length(txl),'xticklabels',txl,'yticklabels',txl,'Ydir','reverse'); xtickangle(45);
    set(gca,'xtick',1:length(txl),'ytick',1:length(txl),'xticklabels',txl,'yticklabels',txl,'Ydir','reverse'); xtickangle(45);
%     text(-0.5,sz+1.1,'Average Overlap Index (N = 5 animals)','FontSize',5); 
    set(gca,'Ydir','normal');ytickangle(25);
    box on
    changePosition(gca,[0.0 0 -0.04 0]);
    hc = putColorBar(gca,[0.08 0.07 -0.11 -0.15],{sprintf('%.1f',minI),sprintf('%.1f',maxI)},6,'eastoutside',[0.1 0.08 0.11 0.15]);
    colormap jet
    save_pdf(hf,mData.pdf_folder,sprintf('OI_Map_%d_mean1_conjunctionP.pdf',ntrials),600);
    break;
end


%% Overlap Indices ImageSC
while 1
    sii = [Ar_t_D ArL_t_D Ars_t_D Ar_t_D ArL_t_D Ars_t_D Ar_t_D ArL_t_D Ars_t_D Ar_i_T ArL_i_T Ars_i_T Ar_i_T ArL_i_T Ars_i_T Ar_i_T ArL_i_T Ars_i_T];
    resp = [FDgFR_D_g_T FDgFR_T_g_D FDgFR_Comp FTgFR_D_g_T FTgFR_T_g_D FTgFR_Comp];

%     [OI,mOI,semOI,OI_mat,p_vals,h_vals] = get_overlap_index(resp,0.5,0.05);
    [OIo,mOI,semOI,OI_mato,p_vals,h_vals,all_CI,mCI,semCI,all_CI_mat,uni] = get_overlap_index(resp,0.5,0.05);
    mOI = mCI;
    sz = size(mOI,1);
%     mOI = OI_mat(:,:,4);
    oM = ones(size(mOI));
    mask1 = (triu(oM,0) & tril(oM,0)); mOI(mask1==1) = NaN;
    maxI = max([mOI(:)]);%semOI(:)]);    
    minI = min([mOI(:)]);%semOI(:)]);
    
    mask = tril(NaN(size(mOI)),0); mask(mask==0) = 1;
    ptxl = {'D','D','D','T','T','T','C','C','C','D','D','D','T','T','T','C','C','C'};
    for ii = 1:length(ptxl)
        txl{ii} = sprintf('%s-%s',ptxl{ii},rasterNamesTxt{sii(ii)});
    end

    
%     mOI = mOI .* mask;
    imAlpha=ones(size(mOI));    %imAlpha(isnan(mask))=0.25; 
    imAlpha(mask1 == 1) = 0;
%     ff = makeFigureRowsCols(2020,[10 4 6 1.5],'RowsCols',[1 2],'spaceRowsCols',[0.1 0.01],'rightUpShifts',[0.05 0.13],'widthHeightAdjustment',[-240 -150]);
    hf = get_figure(5,[8 3 3.5 3.5]);
    %
%     axes(ff.h_axes(1));
    im1 = imagesc(mOI,[minI,maxI]);    im1.AlphaData = imAlpha;
    set(gca,'color',0.5*[1 1 1]);    colormap parula;    %axis equal
    format_axes(gca);
    for rr = 1:size(mCI,1)
        for cc = 1:size(mCI,1)
            if rr == cc
                text(rr-0.25,cc,sprintf('%.0f',mCI(rr,cc)),'FontSize',5);
            end
        end
    end
    axis equal
    set_axes_limits(gca,[0.5 sz+0.5],[0.5 sz+0.5]);
%     set(gca,'xtick',1:length(txl),'ytick',1:length(txl),'xticklabels',txl,'yticklabels',txl,'Ydir','reverse'); xtickangle(45);
    set(gca,'xtick',1:length(txl),'ytick',1:length(txl),'xticklabels',txl,'yticklabels',txl,'Ydir','reverse'); xtickangle(45);
%     text(-0.5,sz+1.1,'Average Overlap Index (N = 5 animals)','FontSize',5); 
    set(gca,'Ydir','normal');ytickangle(25);
    box on
    changePosition(gca,[0.0 0 -0.04 0]);
    hc = putColorBar(gca,[0.08 0.07 -0.11 -0.15],{sprintf('%.1f',minI),sprintf('%.1f',maxI)},6,'eastoutside',[0.1 0.08 0.11 0.15]);
    colormap jet
    save_pdf(hf,mData.pdf_folder,sprintf('OI_Map_%d_mean1_conjunctionP.pdf',ntrials),600);
    break;
end


%% agglomerative hierarchical clustering
while 1
    mOI1 = mOI;
%     mOI1(isnan(mOI1)) = 1;
%     Di = pdist(mOI1);
    Di = pdist(mOI1,@naneucdist);
%     tree = linkage(mOI1,'average','euclidean');
    tree = linkage(Di,'average');
    figure(hf);clf
    [H,T,TC] = dendrogram(tree,'Orientation','top','ColorThreshold','default');
    hf = gcf;
    set(hf,'Position',[7 3 6.97 1]);
    set(H,'linewidth',1);
    set(gca,'xticklabels',txl(TC));xtickangle(45);
    format_axes(gca);
    hx = ylabel('Eucledian Distance');changePosition(hx,[0 -0.1 0]);
    changePosition(gca,[0.03 0.0 0.05 0.05]);
    save_pdf(hf,mData.pdf_folder,sprintf('OI_Map_cluster1.pdf'),600);
    %%
    break;
end
%% I want to explore the trial to trial changes when the animal is running a fixed distance
RsDC = combine_rasters_conditions(RsDt);
plotRasters_simplest(RsDC{an,1},find(resp))

RsTC = combine_rasters_conditions(RsTi);
plotRasters_simplest(RsTC{1,1},[])
%% combine the distance and time rasters horizontally and see individual cell rasters
RsDTC = combine_rasters_horizontally([RsDC,RsTC]);
ccs = cell_list_op(respDT.inh(1,:),[],'or',1);
plotRasters_simplest(RsDTC{an,1},find(ccs{an}));
%%
plotRasters_simplest(RsDTC{an,1},find(respA));
%% find mean over the 30 trials from conditions 3 4 5 and see population vector
mRDTC = calc_mean_rasters(RsDTC,[]); 
%%
an = 1;
Rs = RsDTC(an) ;mR = mRDTC(an);
ccs = cell_list_op(respDTC.resp(an,1),[],'or',1);
ccs =  cell_list_op(ccs,{respA},'or');
% ccs = {respA};
ff = makeFigureRowsCols(107,[1 0.5 4 0.5],'RowsCols',[2 1],...
    'spaceRowsCols',[0 0.01],'rightUpShifts',[0.11 0.13],'widthHeightAdjustment',...
    [-50 -80]);    set(gcf,'color','w');  set(gcf,'Position',[10 3 3.5 5]);
[CRc,aCRc,mRR] = find_population_vector_corr(Rs,mR,ccs,0);
ff = show_population_vector_and_corr_combined(mData,ff,Rs,mRR,CRc,[],[]);
axes(ff.h_axes(2,1));
set(gca,'xtick',[1 49 51 100],'xticklabels',{'0','','','15'});
text(15,-4,'Distance (cm)','FontSize',6); text(65,-4,'Time (s)','FontSize',6);
text(44,-4,'150','FontSize',6); text(50,-4,'0','FontSize',6);
%%
an = 1;
rasters = RsDC{an}.sp_rasters1;
Dur = RsDC{an}.duration1;
rastersT = RsTC{an}.sp_rasters1;
DurT = RsTC{an}.duration1;
nbins = 4;
nShuffle = 500;
numCells = size(rasters,3);
parfor ii = 1:numCells
    rng(3,'twister');
    [OD(ii),~] = info_metrics_S_onlyMI(rasters(:,:,ii),[],nbins,Dur,nShuffle);
    rng(3,'twister');
    [OT(ii),~] = info_metrics_S_onlyMI(rastersT(:,:,ii),[],nbins,DurT,nShuffle);
end
%%
for ii = 1:numCells
    zMIsD(ii) = OD(ii).ShannonMI_Zsh;
    zMIsT(ii) = OT(ii).ShannonMI_Zsh;
end
%%
an = 1;
rasters = RsDC{an}.sp_rasters1;
rastersT = RsTC{an}.sp_rasters1;
rastersDT = RsDTC{an}.sp_rasters1;
numCells = size(rasters,3);
p = zeros(numCells,1);resp = logical(p);
pT = p; respT = logical(p); pDT = p; respDT = logical(p);
parfor ii = 1:numCells
    p(ii) = anova1(rasters(:,:,ii),1:size(rasters(:,:,ii),2),'off');    
    if p(ii) < 0.05
        resp(ii) = 1;
    end
    pT(ii) = anova1(rastersT(:,:,ii),1:size(rastersT(:,:,ii),2),'off');
    if pT(ii) < 0.05
        respT(ii) = 1;
    end
    pDT(ii) = anova1(rastersDT(:,:,ii),1:size(rastersDT(:,:,ii),2),'off');
    if pDT(ii) < 0.05
        respDT(ii) = 1;
    end
    
end
respA = resp | respT | respDT;
% plotRasters_simplest(RsDTC{an,1},find(respA));
%% fitting of 2D gaussian on the raster plots
statsetfitnlm = statset('fitnlm');
statsetfitnlm.MaxIter = 1000;
statsetfitnlm.TolFun = 1e-10;
% statsetfitnlm.Display = 'iter';
statsetfitnlm.TolX = statsetfitnlm.TolFun;
statsetfitnlm.UseParallel = 1;
% statsetfitnlm.RobustWgtFun = 'welsch';

%%
for ii = 1:numCells
    if ~respDT(ii)
        continue;
    end
%     raster = double(rasters(:,:,ii) > 0);
    raster = rastersDT(:,:,ii);
    rasterFilt = imgaussfilt(raster,[1 2]);
    [rasterF,mdl,coeff_rs] = do_gauss_fit2D(rasterFilt,statsetfitnlm,[1 1]);
    figure(10000);clf;
    subplot 131; imagesc(raster);set(gca,'YDir','normal')
    subplot 132; imagesc(rasterFilt);set(gca,'YDir','normal')
    subplot 133; imagesc(rasterF);set(gca,'YDir','normal')
    title(ii);
    pause(1);
%     pT(ii) = anova1(rastersT(:,:,ii),1:size(rastersT(:,:,ii),2),'off');    
%     pDT(ii) = anova1(rastersDT(:,:,ii),1:size(rastersDT(:,:,ii),2),'off');    
end
%%
respA = resp | respT | respDT;
plotRasters_simplest(RsDTC{an,1},find(respA));


