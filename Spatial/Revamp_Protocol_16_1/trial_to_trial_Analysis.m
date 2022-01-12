function trial_to_trial_Analysis

%% find spatial trial to trial correlation
while 1
    si = [Ar_t_D Ar_i_T ArT_t_D ArT_i_T Ars_t_D Ars_i_T ArL_t_D ArL_i_T Arss_t_D Arss_i_T ArC_t_D ArC_i_T Arsss_t_D Arsss_i_T Tb_t_T Lb_t_T Ab_t_T];
%     si = [Ar_i_T ArL_i_T Ars_i_T];
    Rs = o.Rs(:,si);mR = o.mR(:,si);
    for cn = 1:length(si)
    trials = mat2cell([1:10]',ones(size([1:10]')));
    RsC = repmat(Rs(:,cn),1,10);
    mRsCT = cell(size(RsC,1),length(trials));
    for ii = 1:length(trials)
        ii;
        [mRsCT(:,ii),~] = calc_mean_rasters(RsC(:,1),trials{ii});
    end
    allmRsT{cn} = mRsCT;
    end
    disp('Done');
    %%
    break;
end

%% Overlap Indices ImageSC spatial
while 1
    allresp = []; ind = 1;
    for cn = 1:length(si)
        mRsCT = allmRsT{cn};
        for rr = 1:size(mRsCT,1)
            for cc = 1:size(mRsCT,2)
                this_mat = mRsCT{rr,cc};
                resp{rr,cc} = sum(this_mat,2) > 0;
                if rr == 1
                    txl{ind} = sprintf('C%dT%d',cn,cc);
                    ind = ind + 1;
                end
            end
        end
        allresp = [allresp resp];
    end
    an  = 1:3;
    %%
    [OI,mOI,semOI,OI_mat,p_vals,h_vals] = get_overlap_index(allresp(an,:),0.5,0.05);
    %%
%     mOI = OI{3};
    sz = size(mOI,1);
%     mOI = OI_mat(:,:,4);
    oM = ones(size(mOI));
    mask1 = (triu(oM,0) & tril(oM,0)); mOI(mask1==1) = NaN;
    maxI = max([mOI(:);semOI(:)]);    
    minI = min([mOI(:);semOI(:)]);
    
    mask = tril(NaN(size(mOI)),0); mask(mask==0) = 1; 
%     mOI = mOI .* mask;
    imAlpha=ones(size(mOI));    %imAlpha(isnan(mask))=0.25; 
    imAlpha(mask1 == 1) = 0;
%     ff = makeFigureRowsCols(2020,[10 4 6 1.5],'RowsCols',[1 2],'spaceRowsCols',[0.1 0.01],'rightUpShifts',[0.05 0.13],'widthHeightAdjustment',[-240 -150]);
    hf = get_figure(6,[8 1 7 7]);
    %
%     axes(ff.h_axes(1));
    im1 = imagesc(mOI,[minI,maxI]);    im1.AlphaData = imAlpha;
    for ii = 1:16
        plot([(10.5+((ii-1)*10)) (10.5+((ii-1)*10))],[0 170.5],'w','linewidth',0.1); 
        plot([0 170.5],[(10.5+((ii-1)*10)) (10.5+((ii-1)*10))],'w','linewidth',0.1); 
    end
    
%     plot([10.5 10.5],[0 30.5],'r'); plot([20.5 20.5],[0 30.5],'r');
%     plot([0 30.5],[10.5 10.5],'r'); plot([0 30.5],[20.5 20.5],'r');
    set(gca,'color',0.5*[1 1 1]);    colormap jet;    %axis equal
    format_axes(gca);
    set_axes_limits(gca,[0.5 sz+0.5],[0.5 sz+0.5]);
%     set(gca,'xtick',1:length(txl),'ytick',1:length(txl),'xticklabels',txl,'yticklabels',txl,'Ydir','reverse'); xtickangle(45);
    text(-0.3,sz+1.1,'Average Overlap Index (N = 5 animals)','FontSize',5); set(gca,'Ydir','normal');%ytickangle(20);
    box on
    changePosition(gca,[0.03 0 -0.04 0]);
    hc = putColorBar(gca,[0.09 0.07 -0.11 -0.15],{sprintf('%d',minI),sprintf('%.1f',maxI)},6,'eastoutside',[0.1 0.05 0.06 0.05]);
%     save_pdf(hf,mData.pdf_folder,sprintf('OI_Map_mean_tu_spatial.pdf'),600);
    %%
    break;
end

%% spatial agglomerative hierarchical clustering
while 1
    mOI1 = mOI;
    mOI1(isnan(mOI1)) = 1;
    Di = pdist(mOI1);
    tree = linkage(Di);
    figure(hf);clf
    [H,T,TC] = dendrogram(tree,'Orientation','right','ColorThreshold','default');
    hf = gcf;
    set(hf,'Position',[7 3 6 6]);
    set(H,'linewidth',1);
    set(gca,'yticklabels',txl(TC));ytickangle(30);
    format_axes(gca);
    hx = xlabel('Eucledian Distance');%changePosition(hx,[-0.051 0 0]);
    changePosition(gca,[0 0.0 -0.05 0.05]);
    save_pdf(hf,mData.pdf_folder,sprintf('OI_Map_cluster_tu_spatial.pdf'),600);
    %%
    break;
end


