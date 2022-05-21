function trial_to_trial_Analysis

%% load data
while 1
    si = [Ar_t_D ArL_t_D Ars_t_D];
    si = [Ar_t_D Ar_i_T ArL_t_D ArL_i_T Ars_t_D Ars_i_T];
    si = [Lb_T ArL_L_T Lbs_T Ab_T Abs_T Ar_t_D ArL_t_D Ars_t_D Ar_i_T ArL_i_T Ars_i_T];
%     si = [Lb_T ArL_L_T Lbs_T Ab_t_T Ab_i_T Abs_t_T Abs_i_T Ar_t_D Ar_i_T ArL_t_D ArL_i_T Ars_t_D Ars_i_T];
    si = [Ar_t_D Ar_i_T ArL_t_D ArL_i_T Ars_t_D Ars_i_T];
%     si = [Ar_t_D Ar_i_T ArL_t_D ArL_i_T Ars_t_D Ars_i_T Ar_t_T Ar_i_D ArL_t_T ArL_i_D Ars_t_T Ars_i_D];
    Rs = o.Rs(:,si);mR = o.mR(:,si); RsG = Rs; siG = si;
    avgProps = get_props_Rs(Rs,[10,100]); respM = avgProps.good_FR;
    for cn = 1:length(si)
        trials = mat2cell([1:10]',ones(size([1:10]')));
        RsC = repmat(Rs(:,cn),1,10);
        mRsCT = cell(size(RsC,1),length(trials));
        for ii = 1:length(trials)
            ii;
            [mRsCT(:,ii),~] = calc_mean_rasters(RsC(:,1),trials{ii});
        end
        allmRsT{cn} = mRsCT;
        allRsC{cn} = RsC;
    end
    lnsi = length(si);
    disp('Done');
    %%
    break;
end

%% find remapping (code will take long time to run)
while 1
    si = [Ar_t_D ArL_t_D Ars_t_D];
    si = [Ar_t_D Ar_i_T ArL_t_D ArL_i_T Ars_t_D Ars_i_T];
    si = [Lb_T ArL_L_T Lbs_T Ab_T Abs_T Ar_t_D ArL_t_D Ars_t_D Ar_i_T ArL_i_T Ars_i_T];
%     si = [Lb_T Ab_T Ar_t_D Ar_i_T ArL_t_D ArL_i_T Ars_t_D Ars_i_T Lbs_T Abs_T];
%     si = [Ar_i_T ArL_i_T Ars_i_T];
    Rs = o.Rs(:,si);mR = o.mR(:,si); RsG = Rs;
    avgProps = get_props_Rs(Rs,[0,40]); respG{1} = avgProps.good_FR; avgProps = get_props_Rs(Rs,[40,70]); respG{2} = avgProps.good_FR;
    avgProps = get_props_Rs(Rs,[70,100]); respG{3} = avgProps.good_FR; avgProps = get_props_Rs(Rs,[40,100]); respG{4} = avgProps.good_FR;
    for cn = 7%:length(si)
        trials = mat2cell([1:10]',ones(size([1:10]')));
        RsC = repmat(Rs(:,cn),1,10);
        RsC = find_responsive_rasters(RsC,1:10);
        mRsCT = allmRsT{cn};
        for pn = 1:4
            [cn pn]
            out_C{cn,pn} = find_population_vector_corr_remap(RsC,mRsCT,respG{pn});
        end
    end
    %%
    break;
end

%% average correlation of all animals
while 1
    trialsRemap = out_C{7,1};
    meanPVcorr = trialsRemap.mean_PV_corr;
    bigMeanPVCorr = [];
    for rr = 1:size(meanPVcorr,1)
        tempBM = [];
        for cc = 1:size(meanPVcorr,2)
            tempBM = [tempBM meanPVcorr{rr,cc}];
        end
        bigMeanPVCorr = [bigMeanPVCorr;tempBM];
    end
    hf = get_figure(6,[8 3 7 7]);
    imagesc(bigMeanPVCorr);
    %
    %     axes(ff.h_axes(1));
    plot([10.5 10.5],[0 30.5],'r'); plot([20.5 20.5],[0 30.5],'r');
    plot([0 30.5],[10.5 10.5],'r'); plot([0 30.5],[20.5 20.5],'r');
    set(gca,'color',0.5*[1 1 1]);    colormap parula;    %axis equal
    format_axes(gca);
    set_axes_limits(gca,[0.5 sz+0.5],[0.5 sz+0.5]);
    set(gca,'xtick',1:length(txl),'ytick',1:length(txl),'xticklabels',txl,'yticklabels',txl,'Ydir','reverse'); xtickangle(45);
    text(-0.3,sz+1.1,'Average Overlap Index (N = 5 animals)','FontSize',5); set(gca,'Ydir','normal');ytickangle(20);
    box on
    changePosition(gca,[0.03 0 -0.04 0]);
    hc = putColorBar(gca,[0.09 0.07 -0.11 -0.15],{sprintf('%d',minI),sprintf('%.1f',maxI)},6,'eastoutside',[0.1 0.05 0.06 0.05]);
    
%     ff = makeFigureRowsCols(106,[1 0.5 9 9],'RowsCols',[10 10],...
%         'spaceRowsCols',[0.02 0.02],'rightUpShifts',[0.1 0.1],'widthHeightAdjustment',...
%         [-40 -40]);
%     set(gcf,'color','w');
%     set(gcf,'Position',[5 1 9 9]);
%     ff = show_remapping_corr_plots(mData,ff,trialsRemap.mean_PV_corr,trialsRemap.xs,[]);
%     save_pdf(ff.hf,mData.pdf_folder,sprintf('remap_corr_C.pdf'),600);
    %%
    break;
end

%% Overlap Indices ImageSC all
while 1
    % for sensory air responses
    siSeA = [Ab_T Abs_T];
    Rs = o.Rs(:,siSeA);mR = o.mR(:,siSeA);
    ntrials = 50;
    pSeA = get_props_Rs(Rs,ntrials);
    respSeA{1} = pSeA.good_FR_and_exc;     respSeA{2} = pSeA.good_FR_and_inh;    respSeA{3} = pSeA.good_FR_and_untuned;
    %
    siSeL = [Lb_T ArL_L_T Lbs_T];
    Rs = o.Rs(:,siSeL);mR = o.mR(:,siSeL);
    ntrials = 50;
    pSeL = get_props_Rs(Rs,ntrials);
    respSeL{1} = pSeL.good_FR_and_exc;     respSeL{2} = pSeL.good_FR_and_inh;    respSeL{3} = pSeL.good_FR_and_untuned;
    
    avgProps = get_props_Rs(RsG,[50,100]); respG = avgProps.good_FR;
    an  = 1:5; eic = 1; sp = 0; intersect_with_global = 0;
    
    allresp = []; ind = 1;
    all_peakL = [];
    for cn = 1:length(si)
        mRsCT = allmRsT{cn};
        resp = []; peak_locations = [];
        for rr = 1:size(mRsCT,1)
            for cc = 1:size(mRsCT,2)
                this_mat = mRsCT{rr,cc};
                [~,peakL] = max(this_mat,[],2);
%                 size_tmat(rr,cc) = size(this_mat,2);
                resp{rr,cc} = sum(this_mat,2) > 0;
                if intersect_with_global
                    resp{rr,cc} = resp{rr,cc} & respG{rr,cn};
                end
                if sp == 1
                    if cn == 1
                        respSe = respSeL{eic}; resp{rr,cc} = resp{rr,cc} & respSe{rr,1};
                    end
                    if cn == 2
                        respSe = respSeL{eic}; resp{rr,cc} = resp{rr,cc} & respSe{rr,2};
                    end
                    if cn == 3
                        respSe = respSeL{eic}; resp{rr,cc} = resp{rr,cc} & respSe{rr,3};
                    end
                    if cn == 4
                        respSe = respSeA{eic}; resp{rr,cc} = resp{rr,cc} & respSe{rr,1};
                    end
                    if cn == 5
                        respSe = respSeA{eic}; resp{rr,cc} = resp{rr,cc} & respSe{rr,2};
                    end
                end
                peakL(~resp{rr,cc}) = NaN;
                peak_locations{rr,cc} = peakL;
                if rr == 1
                    txl{ind} = sprintf('C%dT%d',cn,cc);
                    ind = ind + 1;
                end
            end
        end
        allresp = [allresp resp]; all_peakL = [all_peakL peak_locations];
    end
    i_allresp = cell_list_op(allresp,[],'not');
%     [OI,mOI,semOI,OI_mat,p_vals,h_vals] = get_overlap_index(allresp(an,:),0.5,0.05,0);
%     [AS,mAS,semAS,AS_mat] = get_average_shift(all_peakL(an,:));
%     [OI,mOI,semOI,OI_mat,mOIM,semOIM,OIM_mat] = get_average_shift(all_peakL(an,:));
%     [OI,mOI,semOI,OI_mat,mOIM,semOIM,OIM_mat] = get_resp_to_non_resp(all_peakL(an,:));
    [OIo,mOI,semOI,OI_mato,p_vals,h_vals,all_CI,mCI,semCI,all_CI_mat,uni] = get_overlap_index(allresp,0.5,0.05);
    mOI = mCI; semOI = semCI;
    break;
end
%% heatmap of conjunction all trials
while 1
%     mOI = mOIM; semOI = semOIM;
%     [OI,mOI,semOI,OI_mat,p_vals,h_vals] = get_overlap_index(i_allresp(an,:),0.5,0.05);
    sz = size(mOI,1);
%     mOI = OI_mat(:,:,4);
    oM = ones(size(mOI));
    mask1 = (triu(oM,0) & tril(oM,0)); mOI(mask1==1) = NaN;
    maxI = max([mOI(:);semOI(:)]);
    minI = min([mOI(:);semOI(:)]);
    
%     mask = tril(NaN(size(mOI)),0); mask(mask==0) = 1; 
%     mOI = mOI .* mask;
    imAlpha=ones(size(mOI));    %imAlpha(isnan(mask))=0.25; 
    imAlpha(mask1 == 1) = 0;
%     ff = makeFigureRowsCols(2020,[10 4 6 1.5],'RowsCols',[1 2],'spaceRowsCols',[0.1 0.01],'rightUpShifts',[0.05 0.13],'widthHeightAdjustment',[-240 -150]);
    hf = get_figure(6,[8 5 3.5 3.5]);
    %
%     axes(ff.h_axes(1));
    im1 = imagesc(mOI,[minI,maxI]);    im1.AlphaData = imAlpha;
    plot([10.5 10.5],[0 30.5],'r'); plot([20.5 20.5],[0 30.5],'r');
    plot([0 30.5],[10.5 10.5],'r'); plot([0 30.5],[20.5 20.5],'r');
    set(gca,'color',0.5*[1 1 1]);    colormap parula;    %axis equal
    format_axes(gca);
    set_axes_limits(gca,[0.5 sz+0.5],[0.5 sz+0.5]);
    ttxl = rasterNamesTxt(siG);
    set(gca,'xtick',1:length(ttxl),'ytick',1:length(ttxl),'xticklabels',ttxl,'yticklabels',ttxl,'Ydir','reverse'); xtickangle(45);
    text(-0.3,sz+1.1,'Average Overlap Index (N = 5 animals)','FontSize',5); set(gca,'Ydir','normal');ytickangle(20);
    box on
    changePosition(gca,[0.03 0 -0.04 0]);
    hc = putColorBar(gca,[0.09 0.07 -0.11 -0.15],{sprintf('%.1f',minI),sprintf('%.1f',maxI)},6,'eastoutside',[0.1 0.05 0.06 0.05]);
    colormap jet
    save_pdf(hf,mData.pdf_folder,sprintf('OI_Map_mean_tu_spatial.pdf'),600);
    %%
    break;
end

%% heatmap of conjunction averaged over trial-pairs
while 1
    mOI = [];
    st = 1:10:60; se = 10:10:60;
    for rr = 1:length(st)
        for cc = 1:length(st)
            tempM = mCI(st(rr):se(rr),st(cc):se(cc));
            oM = ones(size(tempM));
            mask1 = (triu(oM,0) & tril(oM,0)); tempM(mask1==1) = NaN;
            mOI(rr,cc) = nanmean(tempM,'A');
        end
    end
    wind = ones(6,6);
    mOI = conv2(wind,tmOI,'same')/100;
%     mOI = mOIM; semOI = semOIM;
%     [OI,mOI,semOI,OI_mat,p_vals,h_vals] = get_overlap_index(i_allresp(an,:),0.5,0.05);
    sz = size(mOI,1);
%     mOI = OI_mat(:,:,4);
    oM = ones(size(mOI));
%     mask1 = (triu(oM,0) & tril(oM,0)); mOI(mask1==1) = NaN;
    maxI = max([mOI(:);semOI(:)]);
    minI = min([mOI(:);semOI(:)]);
    
%     mask = tril(NaN(size(mOI)),0); mask(mask==0) = 1; 
%     mOI = mOI .* mask;
    imAlpha=ones(size(mOI));    %imAlpha(isnan(mask))=0.25; 
%     imAlpha(mask1 == 1) = 0;
%     ff = makeFigureRowsCols(2020,[10 4 6 1.5],'RowsCols',[1 2],'spaceRowsCols',[0.1 0.01],'rightUpShifts',[0.05 0.13],'widthHeightAdjustment',[-240 -150]);
    hf = get_figure(6,[8 5 3.5 3.5]);
    %
%     axes(ff.h_axes(1));
    im1 = imagesc(mOI,[minI,maxI]);    im1.AlphaData = imAlpha;
    plot([10.5 10.5],[0 30.5],'r'); plot([20.5 20.5],[0 30.5],'r');
    plot([0 30.5],[10.5 10.5],'r'); plot([0 30.5],[20.5 20.5],'r');
    set(gca,'color',0.5*[1 1 1]);    colormap parula;    %axis equal
    format_axes(gca);
    set_axes_limits(gca,[0.5 sz+0.5],[0.5 sz+0.5]);
    ttxl = rasterNamesTxt(siG);
    set(gca,'xtick',1:length(ttxl),'ytick',1:length(ttxl),'xticklabels',ttxl,'yticklabels',ttxl,'Ydir','reverse'); xtickangle(45);
    text(-0.3,sz+1.1,'Average Overlap Index (N = 5 animals)','FontSize',5); set(gca,'Ydir','normal');ytickangle(20);
    box on
    changePosition(gca,[0.03 0 -0.04 0]);
    hc = putColorBar(gca,[0.09 0.07 -0.11 -0.15],{sprintf('%.1f',minI),sprintf('%.1f',maxI)},6,'eastoutside',[0.1 0.05 0.06 0.05]);
    colormap jet
    save_pdf(hf,mData.pdf_folder,sprintf('OI_Map_mean_tu_spatial.pdf'),600);
    %%
    break;
end

%% spatial agglomerative hierarchical clustering
while 1
    mOI1 = mOI;
    mOI1(isnan(mOI1)) = 1;
    Di = pdist(mOI1);
    tree = linkage(Di,'average');
    figure(hf);clf
    [H,T,TC] = dendrogram(tree,'Orientation','right','ColorThreshold','default');
    hf = gcf;
    set(hf,'Position',[7 5 2 2]);
    set(H,'linewidth',1);
	ttxl = rasterNamesTxt(siG);
    ttxl = ttxl(TC);
    set(gca,'yticklabels',ttxl);ytickangle(30);
    format_axes(gca);
    hx = xlabel('Eucledian Distance');%changePosition(hx,[-0.051 0 0]);
    changePosition(gca,[0 0.0 -0.05 0.05]);
    save_pdf(hf,mData.pdf_folder,sprintf('OI_Map_cluster_tu_spatial.pdf'),600);
    %%
    break;
end

%% rm anova of one row
while 1
%     rowN = 1;
%     for an = 1:size(OIM_mat,3)
%         var(an,:) = OIM_mat(rowN,:,an);
%     end
    var = squeeze(mean(OIM_mat,1))';
    [within,dvn,xlabels] = make_within_table({'Cond'},[13]);
    dataT = make_between_table({var},dvn);
    ra = RMA(dataT,within);
    ra.ranova
    %%
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels,extras] = get_vals_for_bar_graph_RMA(mData,ra,{'Cond','hsd'},[1 1 1]);
%     xdata = xdataG;
    ptab = 0;
    if ptab h(h==1) = 0; end
    hf = get_figure(5,[8 7 1.75 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = colors;
    
    if ptab
    [hbs,maxY] = plot_bars_p_table(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k','ptable',extras.pvalsTable,...
        'BaseValue',0.01,'xdata',xdatag,'sigFontSize',7,'sigAsteriskFontSize',4,'barWidth',0.5);
    else
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'BaseValue',0.01,'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'yspacing',0.1);
    end
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    txl = rasterNamesTxt(si); 
    xticks = xdata; xticklabels = txl;
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 0.1]); xtickangle(45);
    changePosition(gca,[0.02 -0.01 0.06 0.01]); put_axes_labels(gca,{[],[0 0 0]},{'Avg. Correlation',[0 0 0]});
    %%
    break;
end

%%
while 1   %%
    [OIo,mOIo,semOIo,OI_mato,p_vals,h_vals,CI,mCI,semCI,CI_mat] = get_overlap_index(allresp(an,:),0.5,0.05);
%     [OI,mOI,semOI,OI_mat,p_vals,h_vals] = get_overlap_index(allresp,0.5,0.05,0);
%     mOI = OI{4}; semOI = semOI;
    mOI = mCI; semOI = semCI;
%     [OI,mOI,semOI,OI_mat,p_vals,h_vals] = get_overlap_index(i_allresp(an,:),0.5,0.05);
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
    hf = get_figure(6,[8 3 3.5 3.5]);
    hf = get_figure(6,[8 3 7 7]);
    %
%     axes(ff.h_axes(1));
    im1 = imagesc(mOI,[minI,maxI]);    im1.AlphaData = imAlpha;
    
    for ii = 1:12
        plot([(10.5+((ii-1)*10)) (10.5+((ii-1)*10))],[0 130.5],'w','linewidth',0.1); 
        plot([0 130.5],[(10.5+((ii-1)*10)) (10.5+((ii-1)*10))],'w','linewidth',0.1); 
    end
    
%     plot([10.5 10.5],[0 130.5],'r'); plot([30.5 30.5],[0 130.5],'r'); plot([90.5 90.5],[0 130.5],'r'); plot([100.5 100.5],[0 130.5],'r'); plot([120.5 120.5],[0 130.5],'r');

%     plot([0 130.5],[10.5 10.5],'r'); plot([0 130.5],[30.5 30.5],'r'); plot([0 130.5],[90.5 90.5],'r'); plot([0 130.5],[100.5 100.5],'r'); plot([0 130.5],[120.5 120.5],'r');
    
%     plot([90.5 90.5],[0 130.5],'r'); plot([90.5 90.5],[0 130.5],'r');
%     plot([0 130.5],[10.5 10.5],'r'); plot([0 130.5],[30.5 30.5],'r');

    set(gca,'color',0.5*[1 1 1]);    colormap parula;    %axis equal
    format_axes(gca);
    set_axes_limits(gca,[0.5 sz+0.5],[0.5 sz+0.5]);
    ttxl = rasterNamesTxt(siG);
    xtickvals = 5:10:130;%[5 15 25 60 100 115 125];
%     si = [Lb_T Ab_t_T Ab_i_T Ar_t_D Ar_i_T ArL_t_D ArL_i_T Ars_t_D Ars_i_T Lbs_T Abs_t_T Abs_i_T ArL_L_T];
    xticklabels = {'Lb','Ab-t','Ab-i','Ar-t','Ar-i','ArL-t','ArL-i','Ar*-t','Ar*-i','Lb*','Ab*-t','Ab*-i','ArL-L'};

    set(gca,'xtick',xtickvals,'ytick',xtickvals,'xticklabels',xticklabels,'yticklabels',xticklabels,'Ydir','normal'); xtickangle(45);
%     text(-0.3,sz+1.1,'Average Overlap Index (N = 5 animals)','FontSize',5); set(gca,'Ydir','normal');ytickangle(20);
    box on
    changePosition(gca,[0.03 0 -0.04 0]);
    hc = putColorBar(gca,[0.09 0.07 -0.11 -0.15],{sprintf('%d',minI),sprintf('%.1f',maxI)},6,'eastoutside',[0.1 0.05 0.06 0.05]);
    colormap jet
    save_pdf(hf,mData.pdf_folder,sprintf('OI_Map_mean_tu_spatial.pdf'),600);
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
    set(hf,'Position',[7 5 2 2]);
    set(H,'linewidth',1);
	ttxl = rasterNamesTxt(siG);
    ttxl = ttxl(TC);
    set(gca,'yticklabels',ttxl);ytickangle(30);
    format_axes(gca);
    hx = xlabel('Eucledian Distance');%changePosition(hx,[-0.051 0 0]);
    changePosition(gca,[0 0.0 -0.05 0.05]);
    save_pdf(hf,mData.pdf_folder,sprintf('OI_Map_cluster_tu_spatial.pdf'),600);
    %%
    break;
end

%% may be useful but telling the same story
while 1
    si = [Lb_T ArL_L_T Lbs_T Ab_T Abs_T Ar_t_D ArL_t_D Ars_t_D Ar_i_T ArL_i_T Ars_i_T];
%     si = [Lb_T Ab_T Ar_t_D Ar_i_T ArL_t_D ArL_i_T Ars_t_D Ars_i_T Lbs_T Abs_T];
%     si = [Ar_i_T ArL_i_T Ars_i_T];
    Rs = o.Rs(:,si);mR = o.mR(:,si); RsG = Rs;
    avgProps = get_props_Rs(Rs,[50,100]); respM = avgProps.good_FR;
    respCL = [];
    corrCond = [];
    meanCorrCond = [];
    for rr = 1:size(respM,1)
        tempCL = [];
        for cc = 1:size(respM,2)
            tempCL = [tempCL respM{rr,cc}];
        end
        respCL{rr} = tempCL;
        [~,corrCond{rr},cellnums] = findPopulationVectorPlot(tempCL,[]); % order RV1 according to peak firing
        meanCorrCond(:,:,rr) = corrCond{rr};
        
    end
    hf = get_figure(6,[8 3 7 7]);
    meanCorrCond = mean(meanCorrCond,3);
    imagesc(meanCorrCond);
    
    %%
    figure(100);clf;
    imagesc(respCL{5});
    
    
    %%
    break;
end

%% 1 off diagnoal (uniqe between adjacent trials)
while 1
    xticklabels = rasterNamesTxt(si);
    for an = 1:5
        respRV(an,:) = diag(all_CI_mat(:,:,an));
        conjV(an,:) = diag(all_CI_mat(:,:,an),1);
        comp1V(an,:) = diag(uni(:,:,an),1);
        comp2V(an,:) = diag(uni(:,:,an),-1);
    end
    respV = respRV;
    mrespV = mean(respV); semrespV = std(respV)./sqrt(5);
    respTW = reshape(mrespV,10,lnsi); mrespAct = respTW'; 
    respTW = reshape(semrespV,10,lnsi); semrespAct = respTW';
    
    respV = conjV;
    mrespV = mean(respV); semrespV = std(respV)./sqrt(5);
    mrespV(10:10:((lnsi*10)-1)) = NaN; mrespV(60) = NaN; mrespV(isnan(mrespV)) = []; 
    semrespV(10:10:((lnsi*10)-1)) = NaN; semrespV(60) = NaN; semrespV(isnan(semrespV)) = []; 
    respTW = reshape(mrespV,9,lnsi); mconjAct = respTW'; 
    respTW = reshape(semrespV,9,lnsi); semconjAct = respTW';
    
    respV = comp1V;
    mrespV = mean(respV); semrespV = std(respV)./sqrt(5);
    mrespV(10:10:((lnsi*10)-1)) = NaN; mrespV(60) = NaN; mrespV(isnan(mrespV)) = []; 
    semrespV(10:10:((lnsi*10)-1)) = NaN; semrespV(60) = NaN; semrespV(isnan(semrespV)) = []; 
    respTW = reshape(mrespV,9,lnsi); mcomp1Act = respTW'; 
    respTW = reshape(semrespV,9,lnsi); semcomp1Act = respTW';
    
    respV = comp2V;
    mrespV = mean(respV); semrespV = std(respV)./sqrt(5);
    mrespV(10:10:((lnsi*10)-1)) = NaN; mrespV(60) = NaN; mrespV(isnan(mrespV)) = []; 
    semrespV(10:10:((lnsi*10)-1)) = NaN; semrespV(60) = NaN; semrespV(isnan(semrespV)) = []; 
    respTW = reshape(mrespV,9,lnsi); mcomp2Act = respTW'; 
    respTW = reshape(semrespV,9,lnsi); semcomp2Act = respTW';
    
    mrespActL = NaN;  semrespActL = NaN; 
    xtl = {''}; trialsStr = cellfun(@num2str,trials','un',0);
    xaxL = NaN; xticks = []; xtickL =[];
    for ii = 1:lnsi
        mrespActL = [mrespActL mrespAct(ii,:) NaN]; semrespActL = [semrespActL semrespAct(ii,:) NaN];
        xaxL = [xaxL 1:10 NaN];
        xtl = [xtl trialsStr {''}];
    end
    xax = 1:length(mrespActL); 
    for ii = 1:length(mrespActL)
        if xaxL(ii) == 1 || xaxL(ii) == 9
            xticks = [xticks xax(ii)];
            xtickL = [xtickL xtl(ii)];
        end
    end
    rlcolor = [0.75 0.75 0.75];
    hf = figure(100);clf;
    set(hf,'units','inches','position',[5 5 6.9 1.5]);
%     plot(xax,mrespActL,'color',rlcolor);hold on;
%     plot(xlim,[nanmean(mrespActL) nanmean(mrespActL)],'color',rlcolor);hold on;
    plot(xlim,[nanmean(mconjAct(:)) nanmean(mconjAct(:))],'color','m');hold on;
    iii=1;
    theinds = find(isnan(mrespActL));
    for ii = find(isnan(mrespActL))
        plot([ii ii],[7 25],'b-');
        if iii <= lnsi
            text(ii+2,27,sprintf('%s',xticklabels{iii}),'FontSize',6);
            indsS = (theinds(iii)+1):(theinds(iii+1)-1);
%             shadedErrorBar(indsS,mrespActL(indsS),semrespActL(indsS),{'color',rlcolor},0.5);
            plot(indsS(1:9),mconjAct(iii,:),'m');
            shadedErrorBar(indsS(1:9),mconjAct(iii,:),semconjAct(iii,:),{'color','m'},0.5);
            plot(indsS(1:9),mcomp1Act(iii,:),'c');
            shadedErrorBar(indsS(1:9),mcomp1Act(iii,:),semcomp1Act(iii,:),{'color','c'},0.5);
            plot(indsS(1:9),mcomp2Act(iii,:),'k');
            shadedErrorBar(indsS(1:9),mcomp2Act(iii,:),semcomp2Act(iii,:),{'color','k'},0.5);
            iii=iii+1;
        end
    end
    xlim([0 length(mrespActL)+1]); ylim([5 33]);
    xlabel('Trial-Pairs');ylabel('Cells (%)');box off;
    set(gca,'xtick',xticks,'xticklabel',xtickL);
    legs = {'Conjunctive Cells      ','Complementary Cells 1','Complementary Cells 2',[9.5 0.1 31 0.2]}; 
    putLegendH(gca,legs,{'m','c','k'},'sigR',{[],'anova',[],6});
    format_axes(gca);
    changePosition(gca,[-0.08 0.1 0.17 -0.1]);
    save_pdf(hf,mData.pdf_folder,sprintf('tria_to_trial_unique.pdf'),600);
break;
end

%% 1 off diagnoal (uniqe between adjacent trials)
while 1
    %%
    respV = [];
    for an = 1:5
        respV(an,:) = diag(all_CI_mat(:,:,an));
    end
    
    [within,dvn,xlabels] = make_within_table({'cond','Trials'},[lnsi,10]);
    dataT = make_between_table({respV},dvn);
    ra = RMA(dataT,within);
    ra.ranova
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'cond','hsd'},[1.5 1 1]);
%     h(h==1) = 0;
	xdata = make_xdata([lnsi],[1 1.5]);
    hf = get_figure(5,[8 7 2 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = mData.dcolors;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',7,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.15);
%     set(hbs(1),'FaceColor',tcolors{1},'FaceAlpha',0.75); set(hbs(3),'FaceColor',tcolors{2},'FaceAlpha',0.75);
    maxY = maxY + 5;
    ylims = ylim;
    format_axes(gca);
%     htxt = text(0,maxY+6,sprintf('Any Condition (%d\x00B1%d%%),   All Conditions (%d%%)',round(mra),round(semra),round(mrall)),'FontSize',6);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[5 maxY]); format_axes(gca);
    xticks = xdata; 
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[5 10 20 30]); xtickangle(45)
    changePosition(gca,[0.02 0.01 0.05 0]); put_axes_labels(gca,{[],[0 0 0]},{{'Cells (%)'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('responsivity.pdf'),600);
    
    %%
    mrespV = mean(respV); semrespV = std(respV)./sqrt(5);
    respTW = reshape(mrespV,10,lnsi); mrespAct = respTW';
    respTW = reshape(semrespV,10,lnsi); semrespAct = respTW';
    mrespActL = NaN;  semrespActL = NaN; 
    xtl = {''}; trialsStr = cellfun(@num2str,trials','un',0);
    xaxL = NaN; xticks = []; xtickL =[];
    for ii = 1:lnsi
        mrespActL = [mrespActL mrespAct(ii,:) NaN];
        semrespActL = [semrespActL semrespAct(ii,:) NaN];
        xaxL = [xaxL 1:10 NaN];
        xtl = [xtl trialsStr {''}];
    end
    xax = 1:length(mrespActL); 
    for ii = 1:length(mrespActL)
        if xaxL(ii) == 1 || xaxL(ii) == 10
            xticks = [xticks xax(ii)];
            xtickL = [xtickL xtl(ii)];
        end
    end
    
    hf = figure(100);clf;
    set(hf,'units','inches','position',[5 5 6.9 1]);
    plot(xax,mrespActL,'k');hold on;
    plot(xlim,[nanmean(mrespActL) nanmean(mrespActL)],'k');
    iii=1;
    theinds = find(isnan(mrespActL));
    for ii = find(isnan(mrespActL))
        plot([ii ii],[23 56],'b-');
        if iii <= lnsi
        text(ii+2,55,sprintf('%s',xticklabels{iii}),'FontSize',lnsi);
        indsS = (theinds(iii)+1):(theinds(iii+1)-1);
        shadedErrorBar(indsS,mrespActL(indsS),semrespActL(indsS),{},0.5)
        iii=iii+1;
        end
        
    end
    xlim([0 length(mrespActL)+1]); ylim([20 60]);
    xlabel('Trials');ylabel('Cells (%)');box off;
    set(gca,'xtick',xticks,'xticklabel',xtickL);
%     legs = {'Responsive Cells',[9.5 0.1 34 0.2]}; 
%     putLegendH(gca,legs,{'k'},'sigR',{[],'anova',[],6});
    format_axes(gca);
    changePosition(gca,[-0.08 0.1 0.17 -0.1]);
    save_pdf(hf,mData.pdf_folder,sprintf('tria_to_trial_unique.pdf'),600);
    %%
break;
end
%%