function trial_to_trial_Analysis_brake_vs_nobrake_updated_wo_Arb

%% find spatial trial to trial correlation
while 1
    cpCells = evalin('base','cellposeCells');
    cpCells = cpCells';
    trialNums = [1:10];
   si = [Lb Lbs Ab_On Abs_On Ab_Off Abs_Off Ar_On ArL_On Ars_On Ar_Off ArL_Off Ars_Off ArL_L Ar_L Ars_L];
   si = [Lb Ab_On Ab_Off Ar_On Ar_Off ArL_On ArL_Off Ars_On Ars_Off Lbs Abs_On Abs_Off ArL_L Ar_L Ars_L];
%    si = [Lb Ab_On Ab_Off Ar_On Ar_Off ArL_On ArL_Off  Ars_On Ars_Off Lbs Abs_On Abs_Off];
   event_type = {'1-L','2-AOn','2-AOff','3-AOn','3-AOff','4-AOn','4-AOff','5-AOn','5-AOff','6-L','7-AOn','7-AOff'};
   si = [Lb Ab_On Ab_Off Ab_Offc Ar_On Ar_L Ar_Off Ar_Offc ArL_On ArL_L ArL_Off ArL_Offc Ars_On Ars_L Ars_Off Ars_Offc Lbs Abs_On Abs_Off Abs_Offc];
   event_type = {'1-L','2-AOn','2-AOff','2-Arb','3-AOn','3-A','3-AOff','3-Arb','4-AOn','4-AL','4-AOff','4-Arb','5-AOn','5-A','5-AOff','5-Arb','6-L','7-AOn','7-AOff','7-Arb'};
   
%    si = [Lb Ab_On Ab_Off Ar_On Ar_Off ArL_On ArL_L ArL_Off Ars_On Ars_Off Lbs Abs_On Abs_Off];
%    event_type = {'1-L','2-AOn','2-AOff','3-AOn','3-AOff','4-AOn','4-AL','4-AOff','5-AOn','5-AOff','6-L','7-AOn','7-AOff'};
%    
    Rs = o.Rs(:,si);mR = o.mR(:,si); RsG = Rs; siG = si; propsG = get_props_Rs(RsG,[40,100]); respG = propsG.vals;
    avgProps = get_props_Rs(Rs,[40,100]); respM = avgProps.good_FR;
    for cn = 1:length(si)
        trials = mat2cell([1:10]',ones(size([1:10]')));
        trials = mat2cell([trialNums]',ones(size([trialNums]')));
        RsC = repmat(Rs(:,cn),1,10);
        mRsCT = cell(size(RsC,1),length(trials));
        for ii = 1:length(trials)
            ii;
            [mRsCT(:,ii),~] = calc_mean_rasters(RsC(:,1),trials{ii});
        end
        allmRsT{cn} = mRsCT;
        allRsC{cn} = RsC;
    end
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
    avgProps = get_props_Rs(RsG,[50,100]); 
    respG = avgProps.vals;
    an  = 1:5; eic = 1; sp = 0; intersect_with_global = 0; only_global = 0;
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
                if only_global
                    resp{rr,cc} = respG{rr,cn};
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
%             oc(rr,cn) = find_cells_based_on_cluster(cell2mat(resp(rr,:)));
        end
        allresp = [allresp resp]; all_peakL = [all_peakL peak_locations];
        
    end
    i_allresp = cell_list_op(allresp,[],'not');

    allrespOR = cell_list_op(allresp,[],'or',1);
    allrespAND = cell_list_op(allresp,[],'and',1);
    
    pallrespOR = 100*exec_fun_on_cell_mat(allrespOR,'sum')./exec_fun_on_cell_mat(allrespOR,'length');
    pallrespAND = 100*exec_fun_on_cell_mat(allrespAND,'sum')./exec_fun_on_cell_mat(allrespAND,'length');
    
    [mparOR,semparOR] = findMeanAndStandardError(pallrespOR);
    
    disp('Done');
    %% for Reviewer 3
    allrespORp = find_percent(allrespOR,cpCells);
    [mCs,semCs] = findMeanAndStandardError(allrespORp);
    %%
    [OIo,mOIo,semOIo,OI_mato,p_vals,h_vals,all_CI,mCI,semCI,all_CI_mat,uni] = get_overlap_index(allresp,0.5,0.05);
%     [OIo,mOI,semOI,OI_mato,p_vals,h_vals,all_CI,mCI,semCI,all_CI_mat,uni] = get_overlap_index(all_resp_T,0.5,0.05);
     disp('Done');
    %%
    p_allresp_or = find_percent(cell_list_op(allresp,[],'or',1));
    descriptiveStatistics(p_allresp_or)
    %%
    mUni = nanmean(uni,3); mUni1 = tril(mUni,-1) + tril(mUni,-1)'; mUni2 = triu(mUni,1) + triu(mUni,1)'; mmUni = min(mUni(:)); MmUni = max(mUni(:));
    mOI = mCI; semOI = semCI;
    mSel = mCI;
%     mSel = mUni;
    mSel = mOIo; semOI = semOIo;
    mOI = mSel;
    sz = size(mOI,1);
%     mOI = OI_mat(:,:,4);
    oM = ones(size(mOI));
    mask1 = (triu(oM,0) & tril(oM,0)); mOI(mask1==1) = NaN;
    maxI = max([mOI(:);semOI(:)]);
    minI = min([mOI(:);semOI(:)]);
%     minI = 0; maxI = 0.6;
    mask = tril(NaN(size(mOI)),0); mask(mask==0) = 1; 
%     mOI = mOI .* mask;
    imAlpha=ones(size(mOI));    %imAlpha(isnan(mask))=0.25; 
    imAlpha(mask1 == 1) = 0;
%     ff = makeFigureRowsCols(2020,[10 4 6 1.5],'RowsCols',[1 2],'spaceRowsCols',[0.1 0.01],'rightUpShifts',[0.05 0.18],'widthHeightAdjustment',[-240 -150]);
    hf = get_figure(6,[8 3 3.25 3.25]);
%     hf = get_figure(6,[8 3 4 4]);
    %
%     axes(ff.h_axes(1));
    im1 = imagesc(mOI,[minI,maxI]);    im1.AlphaData = imAlpha;
    
    for ii = 1:20
        plot([(10.5+((ii-1)*10)) (10.5+((ii-1)*10))],[0 200.5],'w','linewidth',0.1); 
        plot([0 200.5],[(10.5+((ii-1)*10)) (10.5+((ii-1)*10))],'w','linewidth',0.1); 
    end
%     for ii = [2 6 9 12]   
%         plot([(10.5+((ii-1)*10)) (10.5+((ii-1)*10))],[0 200.5],'k','linewidth',1); 
%         plot([0 200.5],[(10.5+((ii-1)*10)) (10.5+((ii-1)*10))],'k','linewidth',1); 
%     end
    set(gca,'color',0.5*[1 1 1]);    colormap parula;    %axis equal
    format_axes(gca);
    set_axes_limits(gca,[0.5 sz+0.5],[0.5 sz+0.5]);
    ttxl = rasterNamesTxt(siG);
    xtickvals = 5:10:size(mOI,2);%[5 15 25 60 100 115 125];
    xticklabels = rasterNamesTxt(siG);

    set(gca,'xtick',xtickvals,'ytick',xtickvals,'xticklabels',xticklabels,'yticklabels',xticklabels,'Ydir','normal'); xtickangle(45);%ytickangle(45);
    yyaxis right
    set(gca,'ytick',xtickvals,'yticklabels',yticklabels,'tickdir','out');
    box off
    changePosition(gca,[-0.01 0.00 0.037 0.033]);
    hc = putColorBar(gca,[0.1 -0.08 -0.2 0.03],{sprintf('%.1f',minI),sprintf('%.1f',maxI)},6,'northoutside',[0.07 0.09 0.02 0.09]);
    colormap jet
    save_pdf(hf,mData.pdf_folder,sprintf('OI_Map_mean_tu_spatial.pdf'),600);
    
    figdim = 1;
    hf = get_figure(7,[5 2 figdim+0.1 figdim]);
    im1 = imagesc(semOI,[minI,maxI]);    im1.AlphaData = imAlpha;
    set(gca,'color',0.5*[1 1 1]);    colormap parula;    %axis equal
    format_axes(gca);
    set_axes_limits(gca,[0.5 sz+0.5],[0.5 sz+0.5]);
    set(gca,'xtick',1:length(txl),'ytick',1:length(txl),'xticklabels',[],'yticklabels',[],'Ydir','reverse'); xtickangle(75);
    changePosition(gca,[0.0 0 -0.05 0]);
    set(gca,'Ydir','normal');
    box off
    colormap jet
    save_pdf(hf,mData.pdf_folder,sprintf('OI_Map_sem.pdf'),600);
    
     %%
%      mCI1 = mUni1;
%      mCI1 = mUni2;
%      mCI1 = mCI;
     mask1 = (triu(oM,0) & tril(oM,0)); mCI1(mask1==1) = NaN;
     sts = 1:10:size(mOI,2); ses = 10:10:size(mOI,2);
     clear mCIm mUni1m mUni2m;
     for rr = 1:length(siG)
         for cc = 1:length(siG)
%              mOI(rr,cc) = mean(mCI1(sts(rr):ses(rr),sts(cc):ses(cc)),'All');
             mCIm(rr,cc) = mean(mCI(sts(rr):ses(rr),sts(cc):ses(cc)),'All');
             mUni1m(rr,cc) = mean(mUni1(sts(rr):ses(rr),sts(cc):ses(cc)),'All');
             mUni2m(rr,cc) = mean(mUni2(sts(rr):ses(rr),sts(cc):ses(cc)),'All');
         end
     end
     %%
%      mOI = mCI; semOI = semCI;
    ahc_col_th
    sz = size(mOI,1);
%     mOI = OI_mat(:,:,4);
    oM = ones(size(mOI));
    mask1 = (triu(oM,0) & tril(oM,0)); mOI(mask1==1) = NaN;
    maxI = max([mOI(:);semOI(:)]);
    minI = min([mOI(:);semOI(:)]);
%     minI = 0; maxI = 0.6;
    mask = tril(NaN(size(mOI)),0); mask(mask==0) = 1; 
%     mOI = mOI .* mask;
    imAlpha=ones(size(mOI));    %imAlpha(isnan(mask))=0.25; 
    imAlpha(mask1 == 1) = 0;
%     ff = makeFigureRowsCols(2020,[10 4 6 1.5],'RowsCols',[1 2],'spaceRowsCols',[0.1 0.01],'rightUpShifts',[0.05 0.18],'widthHeightAdjustment',[-240 -150]);
    hf = get_figure(6,[8 3 3.5 3.5]);
%     hf = get_figure(6,[8 3 4 4]);
    %
%     axes(ff.h_axes(1));
    im1 = imagesc(mOI,[minI,maxI]);    im1.AlphaData = imAlpha;
    
    for ii = 1:12
        plot([(10.5+((ii-1)*10)) (10.5+((ii-1)*10))],[0 120.5],'w','linewidth',0.1); 
        plot([0 120.5],[(10.5+((ii-1)*10)) (10.5+((ii-1)*10))],'w','linewidth',0.1); 
    end
%     for ii = [2 6 9 12]   
%         plot([(10.5+((ii-1)*10)) (10.5+((ii-1)*10))],[0 200.5],'k','linewidth',1); 
%         plot([0 200.5],[(10.5+((ii-1)*10)) (10.5+((ii-1)*10))],'k','linewidth',1); 
%     end
    set(gca,'color',0.5*[1 1 1]);    colormap parula;    %axis equal
    format_axes(gca);
    set_axes_limits(gca,[0.5 sz+0.5],[0.5 sz+0.5]);
    ttxl = rasterNamesTxt(siG);
    xtickvals = 5:10:120;%[5 15 25 60 100 115 125];
    xticklabels = rasterNamesTxt(siG);

    set(gca,'xtick',xtickvals,'ytick',xtickvals,'xticklabels',xticklabels,'yticklabels',xticklabels,'Ydir','normal'); xtickangle(45);%ytickangle(45);
    yyaxis right
    set(gca,'ytick',xtickvals,'yticklabels',yticklabels,'tickdir','out');
    box off
    changePosition(gca,[-0.01 0.00 0.037 0.033]);
    hc = putColorBar(gca,[0.1 -0.08 -0.2 0.03],{sprintf('%.1f',minI),sprintf('%.1f',maxI)},6,'northoutside',[0.07 0.09 0.02 0.09]);
    colormap jet
    save_pdf(hf,mData.pdf_folder,sprintf('OI_Map_mean_tu_spatial.pdf'),600);
    disp('Done');
    %%
    ahc_col_th = 0.7;
    mOI1 = mCI;
    mOI1(mask1 == 1) = NaN;
%     mOI1(isnan(mOI1)) = 1;
    Di = pdist(mOI1,@naneucdist);
%     tree = linkage(mOI1,'average','euclidean');
    tree = linkage(Di,'average');
    figure(hf);clf
    [H,T,TC] = dendrogram(tree,'Orientation','top','ColorThreshold',ahc_col_th*max(tree(:,3)));
    hf = gcf;
    txl = event_type;
%     set(hf,'Position',[7 3 3.6 1.9]);
    set(H,'linewidth',1);
%     set(gca,'xticklabels',txl(TC));xtickangle(45);
    format_axes(gca);
    hx = ylabel('Eucledian Distance');changePosition(hx,[0 -0.1 0]);
    changePosition(gca,[-0.05 0.0 0.09 0.05]);
    save_pdf(hf,mData.pdf_folder,sprintf('OI_Map_cluster.pdf'),600);
    %%
    break;
end

%% along diagnol (responsiveness)
while 1
mask = diag(mCI);
hf = figure(200);clf;
set(hf,'units','inches','position',[5 5 3.5 1.5]);
plot(1:200,mask,'m');hold on;
mconj = nanmean(mask);
plot([1 200],[mconj mconj],'k');
xlim([1 200]);
for ii = 1:20
    plot([ii*10 ii*10],ylim,'b-');
%     text(ii*10-5,23.5,sprintf('%s',xticklabels{ii}),'FontSize',6);
end
xticks = 5:10:200;
yticks = [5 10 15 mconj 20];
set(gca,'Xtick',xticks,'XTickLabels',xticklabels);xtickangle(30);
ylabel('Cells (%)');box off;
format_axes(gca);
save_pdf(hf,mData.pdf_folder,sprintf('tria_by_trial_responsive.pdf'),600);
break;
end

%% along diagnol (responsiveness) all on the same trial wise
while 1
mask = diag(mCI);
respTW = reshape(mask,10,20);
respTW = respTW';
hf = figure(100);clf;
set(hf,'units','inches','position',[5 5 3 1.5]);
tcolors = mData.dcolors(1:20);
for ii = 1:20
    plot(1:10,respTW(ii,:),'color',mData.dcolors{ii});hold on;
end
mconj = nanmean(respTW); semconj = std(respTW)./sqrt(5);
h = shadedErrorBar(1:10,mconj,semconj,{'color','k','linewidth',1.5},0.5);
xlim([1 12]); ylim([12,25]);
legs = {xticklabels{:},[10.5 0.1 24 0.3]}; putLegend(gca,legs,tcolors,'sigR',{[],'anova',[],6});
ylabel('Cells (%)');box off;xlabel('Trials');
set(gca,'Xtick',1:10);
format_axes(gca);
changePosition(gca,[-0.03 0.01 0.1 0]);
save_pdf(hf,mData.pdf_folder,sprintf('tria_by_trial_responsive.pdf'),600);
break;
end
%% 1 off diagnoal (conjunctive between adjacent trials)
while 1
mask = diag(mCI,1);
mask(10:10:199) = NaN;
mask(200) = NaN;
mconj = nanmean(mask);
hf = figure(100);clf;
set(hf,'units','inches','position',[5 5 3.5 0.75]);
plot(1:200,mask,'m');hold on;
plot([1 200],[mconj mconj],'k');
xlim([1 200]);
for ii = 1:20
    plot([ii*10 ii*10],[4 11],'b--');
%     text(ii*10-5,23.5,sprintf('%s',xticklabels{ii}),'FontSize',6);
end
xticks = 5:10:200;
set(gca,'Xtick',xticks,'XTickLabels',xticklabels);xtickangle(30);
ylabel('Cells (%)');box off;
format_axes(gca);
save_pdf(hf,mData.pdf_folder,sprintf('tria_to_trial_conj.pdf'),600);
break;
end

%% 1 off diagnoal (conjunctive between adjacent trials) trial-pair wise
while 1
mask = diag(mCI,1);
mask(10:10:199) = NaN;
mask(200) = NaN;
mask(isnan(mask)) = [];
respTW = reshape(mask,9,20);
respTW = respTW';
hf = figure(100);clf;
set(hf,'units','inches','position',[5 5 3 1.5]);
tcolors = mData.dcolors(1:20);
for ii = 1:20
    plot(1:9,respTW(ii,:),'color',mData.dcolors{ii});hold on;
end
mconj = nanmean(respTW); semconj = std(respTW)./sqrt(5);
h = shadedErrorBar(1:9,mconj,semconj,{'color','k','linewidth',1.5},0.5);
xlim([1 11]); ylim([4,11]);
legs = {xticklabels{:},[9.5 0.1 11 0.2]}; putLegend(gca,legs,tcolors,'sigR',{[],'anova',[],6});
box off;
set(gca,'Xtick',1:9,'xticklabels',{'1-2','2-3','3-4','4-5','5-6','6-7','7-8','8-9','9-10'});xtickangle(30)
put_axes_labels(gca,{'Trial-Pairs',[-1 0.5 0]},{{'Cells (%)'},[0.31 0 0]});
format_axes(gca);
changePosition(gca,[-0.03 0.1 0.1 -0.1]);
save_pdf(hf,mData.pdf_folder,sprintf('tria_to_trial_conj.pdf'),600);
break;
end

%% conjunctive trials stat
while 1
    pvals = [];
    for trialN = 1:9
    tvals = [];
    for an = 1:5
        vals = diag(all_CI_mat(:,:,an),1);
        inds = trialN:10:199;
        tvals(an,:) = vals(inds);
    end
    
    [within,dvn,xlabels] = make_within_table({'cond'},[20]);
    dataT = make_between_table({tvals},dvn);
    ra = RMA(dataT,within);
    ra.ranova;
    pvals(trialN) = ra.ranova.pValue_sel(3);
    eta2(trialN) = ra.eta2;
        
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'cond','hsd'},[1.5 1 1]);
%     h(h==1) = 0;
	xdata = make_xdata([20],[1 1.5]);
    hf = get_figure(5,[8 7 1.75 1]);
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
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; 
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 10 20 30]); xtickangle(45)
    changePosition(gca,[0.02 0.01 0.05 0]); put_axes_labels(gca,{[],[0 0 0]},{{'Cells (%)'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('trial_pairs_conj_%d.pdf',trialN),600);
    end
    
    break;
end

%% unique cells
while 1
    mOI = mean(uni,3); semOI = std(uni,[],3)/sqrt(5);
    sz = size(mOI,1);
%     mOI = OI_mat(:,:,4);
    oM = ones(size(mOI));
    mask1 = (triu(oM,0) & tril(oM,0)); mOI(mask1==1) = NaN;
    maxI = max([mOI(:);semOI(:)]);
    minI = min([mOI(:);semOI(:)]);
%     minI = 0; maxI = 0.6;
    mask = tril(NaN(size(mOI)),0); mask(mask==0) = 1; 
%     mOI = mOI .* mask;
    imAlpha=ones(size(mOI));    %imAlpha(isnan(mask))=0.25; 
    imAlpha(mask1 == 1) = 0;
%     ff = makeFigureRowsCols(2020,[10 4 6 1.5],'RowsCols',[1 2],'spaceRowsCols',[0.1 0.01],'rightUpShifts',[0.05 0.18],'widthHeightAdjustment',[-240 -150]);
    hf = get_figure(6,[8 3 3.5 3.5]);
%     hf = get_figure(6,[8 3 4 4]);
    %
%     axes(ff.h_axes(1));
    im1 = imagesc(mOI,[minI,maxI]);    im1.AlphaData = imAlpha;
    
    for ii = 1:12
        plot([(10.5+((ii-1)*10)) (10.5+((ii-1)*10))],[0 200.5],'w','linewidth',0.1); 
        plot([0 200.5],[(10.5+((ii-1)*10)) (10.5+((ii-1)*10))],'w','linewidth',0.1); 
    end
%     for ii = [2 6 9 12]   
%         plot([(10.5+((ii-1)*10)) (10.5+((ii-1)*10))],[0 200.5],'k','linewidth',1); 
%         plot([0 200.5],[(10.5+((ii-1)*10)) (10.5+((ii-1)*10))],'k','linewidth',1); 
%     end
    set(gca,'color',0.5*[1 1 1]);    colormap parula;    %axis equal
    format_axes(gca);
    set_axes_limits(gca,[0.5 sz+0.5],[0.5 sz+0.5]);
    ttxl = rasterNamesTxt(siG);
    xtickvals = 5:10:200;%[5 15 25 60 100 115 125];
    xticklabels = rasterNamesTxt(siG);
%     yticklabels = {'ON','ON','ON','ON','ON','OFF','OFF','ON','ON','ON','OFF','OFF','OFF',};

    set(gca,'xtick',xtickvals,'ytick',xtickvals,'xticklabels',xticklabels,'yticklabels',xticklabels,'Ydir','normal'); xtickangle(45);%ytickangle(45);
    yyaxis right
    set(gca,'ytick',xtickvals,'yticklabels',yticklabels,'tickdir','out');
%     text(-0.3,sz+1.1,'Average Overlap Index (N = 5 animals)','FontSize',5); set(gca,'Ydir','normal');ytickangle(20);
    box off
    changePosition(gca,[-0.01 0.00 0.037 0.033]);
    hc = putColorBar(gca,[0.1 -0.08 -0.2 0.03],{sprintf('%.1f',minI),sprintf('%.1f',maxI)},6,'northoutside',[0.07 0.09 0.02 0.09]);
    colormap jet
    save_pdf(hf,mData.pdf_folder,sprintf('OI_Map_mean_tu_spatial_unique.pdf'),600);
    break;
end

%% 1 off diagnoal (uniqe between adjacent trials)
while 1
mask = diag(mOI,-1);
mask(10:10:199) = NaN;
mask(200) = NaN;
mconj = nanmean(mask);
hf = figure(100);clf;
set(hf,'units','inches','position',[5 5 3.5 0.75]);
plot(1:200,mask,'m');hold on;
plot([1 200],[mconj mconj],'k');
xlim([1 200]);
for ii = 1:20
    plot([ii*10 ii*10],ylim,'b--');
%     text(ii*10-5,23.5,sprintf('%s',xticklabels{ii}),'FontSize',6);
end
xticks = 5:10:200;
set(gca,'Xtick',xticks,'XTickLabels',xticklabels);xtickangle(30);
ylabel('Cells (%)');box off;
format_axes(gca);
save_pdf(hf,mData.pdf_folder,sprintf('tria_to_trial_unique.pdf'),600);
break;
end

%% 1 off diagnoal (unique ) trial-pair wise
while 1

mask = diag(mOI,-1);
mask(10:10:199) = NaN;
mask(200) = NaN;
mask(isnan(mask)) = [];
respTW = reshape(mask,9,20);
respTW = respTW';
hf = figure(100);clf;
set(hf,'units','inches','position',[5 5 3 1.5]);
tcolors = mData.dcolors(1:20);
for ii = 1:20
    plot(1:9,respTW(ii,:),'color',mData.dcolors{ii});hold on;
end
mconj = nanmean(respTW); semconj = std(respTW)./sqrt(5);
h = shadedErrorBar(1:9,mconj,semconj,{'color','k','linewidth',1.5},0.5);
xlim([1 11]); ylim([4,20]);
legs = {xticklabels{:},[9.5 0.1 17 0.2]}; putLegend(gca,legs,tcolors,'sigR',{[],'anova',[],6});
box off;
set(gca,'Xtick',1:9,'xticklabels',{'1-2','2-3','3-4','4-5','5-6','6-7','7-8','8-9','9-10'});xtickangle(30)
put_axes_labels(gca,{'Trial-Pairs',[-1 0.5 0]},{{'Cells (%)'},[0.31 0 0]});
format_axes(gca);
changePosition(gca,[-0.03 0.1 0.1 -0.1]);
save_pdf(hf,mData.pdf_folder,sprintf('tria_to_trial_conj.pdf'),600);
break;
end


%% unique trials stat
while 1
    pvals = [];
    for trialN = 1:9
    tvals = [];
    for an = 1:5
        vals = diag(uni(:,:,an),-1);
        inds = trialN:10:199;
        tvals(an,:) = vals(inds);
    end
    
    [within,dvn,xlabels] = make_within_table({'cond'},[20]);
    dataT = make_between_table({tvals},dvn);
    ra = RMA(dataT,within);
    ra.ranova;
    pvals(trialN) = ra.ranova.pValue_sel(3);
    eta2(trialN) = ra.eta2;
        
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'cond','hsd'},[1.5 1 1]);
%     h(h==1) = 0;
	xdata = make_xdata([20],[1 1.5]);
    hf = get_figure(5,[8 7 1.75 1]);
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
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; 
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 10 20 30]); xtickangle(45)
    changePosition(gca,[0.02 0.01 0.05 0]); put_axes_labels(gca,{[],[0 0 0]},{{'Cells (%)'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('trial_pairs_conj_%d_unique.pdf',trialN),600);
    end
    break;
end

%% 1 off diagnoal (disrupted between adjacent trials)
while 1
mask = diag(mOI,1);
mask(10:10:199) = NaN;
mask(200) = NaN;
mconj = nanmean(mask);

masku = diag(mOI,-1);
masku(10:10:199) = NaN;
masku(200) = NaN;
mconju = nanmean(masku);

hf = figure(100);clf;
set(hf,'units','inches','position',[5 5 3.5 0.75]);
plot(1:200,mask,'m');hold on;
plot(1:200,masku,'c');
plot([1 200],[mconj mconj],'k');
xlim([1 200]);
for ii = 1:20
    plot([ii*10 ii*10],ylim,'b--');
%     text(ii*10-5,23.5,sprintf('%s',xticklabels{ii}),'FontSize',6);
end
xticks = 5:10:200;
set(gca,'Xtick',xticks,'XTickLabels',xticklabels);xtickangle(30);
ylabel('Cells (%)');box off;
format_axes(gca);
save_pdf(hf,mData.pdf_folder,sprintf('tria_to_trial_disrupted.pdf'),600);
break;
end

%% 1 off diagnoal (disrupted ) trial-pair wise
while 1

mask = diag(mOI,1);
mask(10:10:199) = NaN;
mask(200) = NaN;
mask(isnan(mask)) = [];
respTW = reshape(mask,9,20);
respTW = respTW';
hf = figure(100);clf;
set(hf,'units','inches','position',[5 5 3 1.5]);
tcolors = mData.dcolors(1:20);
for ii = 1:20
    plot(1:9,respTW(ii,:),'color',mData.dcolors{ii});hold on;
end
mconj = nanmean(respTW); semconj = std(respTW)./sqrt(5);
h = shadedErrorBar(1:9,mconj,semconj,{'color','k','linewidth',1.5},0.5);
xlim([1 11]); ylim([4,20]);
legs = {xticklabels{:},[9.5 0.1 17 0.2]}; putLegend(gca,legs,tcolors,'sigR',{[],'anova',[],6});
box off;
set(gca,'Xtick',1:9,'xticklabels',{'1-2','2-3','3-4','4-5','5-6','6-7','7-8','8-9','9-10'});xtickangle(30)
put_axes_labels(gca,{'Trial-Pairs',[-1 0.5 0]},{{'Cells (%)'},[0.31 0 0]});
format_axes(gca);
changePosition(gca,[-0.03 0.1 0.1 -0.1]);
save_pdf(hf,mData.pdf_folder,sprintf('tria_to_trial_conj.pdf'),600);
break;
end


%% disrupted trials stat
while 1
    pvals = [];
    for trialN = 1:9
    tvals = [];
    for an = 1:5
        vals = diag(uni(:,:,an),1);
        inds = trialN:10:199;
        tvals(an,:) = vals(inds);
    end
    
    [within,dvn,xlabels] = make_within_table({'cond'},[20]);
    dataT = make_between_table({tvals},dvn);
    ra = RMA(dataT,within);
    ra.ranova;
    pvals(trialN) = ra.ranova.pValue_sel(3);
    eta2(trialN) = ra.eta2;
        
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'cond','hsd'},[1.5 1 1]);
%     h(h==1) = 0;
	xdata = make_xdata([20],[1 1.5]);
    hf = get_figure(5,[8 7 1.75 1]);
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
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; 
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 10 20 30]); xtickangle(45)
    changePosition(gca,[0.02 0.01 0.05 0]); put_axes_labels(gca,{[],[0 0 0]},{{'Cells (%)'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('trial_pairs_conj_%d_disrupted.pdf',trialN),600);
    end
    break;
end


%% 1 off diagnoal (uniqe between adjacent trials) 
% this is the one being used
while 1
    
    for an = 1:5
        respRV(an,:) = diag(all_CI_mat(:,:,an));
        conjV(an,:) = diag(all_CI_mat(:,:,an),1);
        comp1V(an,:) = diag(uni(:,:,an),1);
        comp2V(an,:) = diag(uni(:,:,an),-1);
    end
    respV = respRV;
    mrespV = mean(respV); semrespV = std(respV)./sqrt(5);
    respTW = reshape(mrespV,10,20); mrespAct = respTW'; 
    respTW = reshape(semrespV,10,20); semrespAct = respTW';
    
    respV = conjV;
    mrespV = mean(respV); semrespV = std(respV)./sqrt(5);
    mrespV(10:10:199) = NaN; mrespV(200) = NaN; mrespV(isnan(mrespV)) = []; 
    semrespV(10:10:199) = NaN; semrespV(200) = NaN; semrespV(isnan(semrespV)) = []; 
    respTW = reshape(mrespV,9,20); mconjAct = respTW'; 
    respTW = reshape(semrespV,9,20); semconjAct = respTW';
    
    respV = comp1V;
    mrespV = mean(respV); semrespV = std(respV)./sqrt(5);
    mrespV(10:10:199) = NaN; mrespV(200) = NaN; mrespV(isnan(mrespV)) = []; 
    semrespV(10:10:199) = NaN; semrespV(200) = NaN; semrespV(isnan(semrespV)) = []; 
    respTW = reshape(mrespV,9,20); mcomp1Act = respTW'; 
    respTW = reshape(semrespV,9,20); semcomp1Act = respTW';
    
    respV = comp2V;
    mrespV = mean(respV); semrespV = std(respV)./sqrt(5);
    mrespV(10:10:199) = NaN; mrespV(200) = NaN; mrespV(isnan(mrespV)) = []; 
    semrespV(10:10:199) = NaN; semrespV(200) = NaN; semrespV(isnan(semrespV)) = []; 
    respTW = reshape(mrespV,9,20); mcomp2Act = respTW'; 
    respTW = reshape(semrespV,9,20); semcomp2Act = respTW';
    
    mrespActL = NaN;  semrespActL = NaN; 
    xtl = {''}; trialsStr = cellfun(@num2str,trials','un',0);
    xaxL = NaN; xticks = []; xtickL =[];
    for ii = 1:20
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
        plot([ii ii],[4 21],'b-');
        if iii <= 20
            text(ii+2,21,sprintf('%s',xticklabels{iii}),'FontSize',6);
            indsS = (theinds(iii)+1):(theinds(iii+1)-1);
%             shadedErrorBar(indsS,mrespActL(indsS),semrespActL(indsS),{'color',rlcolor},0.5);
            plot(indsS(1:9),mconjAct(iii,:),'m');
            shadedErrorBar(indsS(1:9),mconjAct(iii,:),semconjAct(iii,:),{'color','m'},0.5);
            plot(indsS(1:9),mcomp1Act(iii,:),'m');
            shadedErrorBar(indsS(1:9),mcomp1Act(iii,:),semcomp1Act(iii,:),{'color','c'},0.5);
            plot(indsS(1:9),mcomp2Act(iii,:),'m');
            shadedErrorBar(indsS(1:9),mcomp2Act(iii,:),semcomp2Act(iii,:),{'color','k'},0.5);
            iii=iii+1;
        end
    end
    xlim([0 length(mrespActL)+1]); ylim([3 27]);
    xlabel('Trial-Pairs');ylabel('Cells (%)');box off;
    set(gca,'xtick',xticks,'xticklabel',xtickL);
    legs = {'Conjunctive Cells      ','Complementary Cells 1','Complementary Cells 2',[9.5 0.1 25 0.2]}; 
    putLegendH(gca,legs,{'m','c','k'},'sigR',{[],'anova',[],6});
    format_axes(gca);
    changePosition(gca,[-0.08 0.1 0.17 -0.1]);
    save_pdf(hf,mData.pdf_folder,sprintf('tria_to_trial_unique.pdf'),600);
break;
end
%% Reviewer 1 2nd last comment
% mR_last_trial_brake = find_percent(cell_list_op(allresp(:,[9]),[],'or',1));
% mR_first_trial_no_brake = find_percent(cell_list_op(allresp(:,[11]),[],'or',1));
% % 
% [h,p,cd,stat] = ttest(mR_last_trial_brake,mR_first_trial_no_brake)

mR_last_trial_brake = find_percent(cell_list_op(allresp(:,[20 30 40]),[],'or',1));
mR_first_trial_no_brake = find_percent(cell_list_op(allresp(:,[41 61 71]),[],'or',1));
% 
[h,p,cd,stat] = ttest(mR_last_trial_brake,mR_first_trial_no_brake)

mR_last_trial_brake = find_percent(cell_list_op(allresp(:,[50 70 80]),[],'or',1));
mR_first_trial_no_brake = find_percent(cell_list_op(allresp(:,[81 101 111]),[],'or',1));
% 
[h,p,cd,stat] = ttest(mR_last_trial_brake,mR_first_trial_no_brake)
% %%
% mR_last_trial_brake = [find_percent(cell_list_op(allresp(:,[19 29 39]),[],'or',1)) find_percent(cell_list_op(allresp(:,[20 30 40]),[],'or',1))];
% mR_first_trial_no_brake = [find_percent(cell_list_op(allresp(:,[51 71 81]),[],'or',1)) find_percent(cell_list_op(allresp(:,[52 72 82]),[],'or',1))];
% % 
% % [h,p,cd,stat] = ttest(mR_last_trial_brake,mR_first_trial_no_brake)
% [within,dvn,xlabels] = make_within_table({'Brake','Trials'},[2,2]);
% dataT = make_between_table({[mR_last_trial_brake mR_first_trial_no_brake]},dvn);
% ra = RMA(dataT,within);
% ra.ranova;
% print_for_manuscript(ra)
% 
% [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Trials','hsd'},[1.5 1 1]);
% %     h(h==1) = 0;
% xdata = make_xdata([2],[1 1.5]);
%   hf = get_figure(5,[8 7 1.25 1.25]);
%   % s = generate_shades(length(bins)-1);
%   tcolors = mData.dcolors;
%   [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
%       'ySpacing',3,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
%       'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.15);
% %     set(hbs(1),'FaceColor',tcolors{1},'FaceAlpha',0.75); set(hbs(3),'FaceColor',tcolors{2},'FaceAlpha',0.75);
%   maxY = maxY + 1;
%   ylims = ylim;
%   format_axes(gca);
% %     htxt = text(0,maxY+6,sprintf('Any Condition (%d\x00B1%d%%),   All Conditions (%d%%)',round(mra),round(semra),round(mrall)),'FontSize',6);
%   set_axes_limits(gca,[0.35 xdata(end)+.65],[5 maxY]); format_axes(gca);
%   xticks = xdata; 
%   set(gca,'xtick',xticks,'xticklabels',{'Trials 9,10','Trials 1,2'},'ytick',[5 15 25]); xtickangle(45)
%   changePosition(gca,[0.0 0.01 0.05 0]); put_axes_labels(gca,{[],[0 0 0]},{{'Estimated','Marginal Means'},[0 0 0]});
%   save_pdf(hf,mData.pdf_folder,sprintf('responsivity.pdf'),600);
%% Reviewer 1 last comment
all_resp_T = allresp(:,1:10);csize = size(all_resp_T,2)+1;
for ttt = 1:10
  all_resp_T(:,csize) = allresp(:,10+ttt);csize = size(all_resp_T,2)+1;
  all_resp_T(:,csize) = allresp(:,20+ttt);csize = size(all_resp_T,2)+1;
  all_resp_T(:,csize) = allresp(:,30+ttt);csize = size(all_resp_T,2)+1;
end
for ttt = 1:10
  all_resp_T(:,csize) = allresp(:,40+ttt);csize = size(all_resp_T,2)+1;
  all_resp_T(:,csize) = allresp(:,50+ttt);csize = size(all_resp_T,2)+1;
  all_resp_T(:,csize) = allresp(:,60+ttt);csize = size(all_resp_T,2)+1;
  all_resp_T(:,csize) = allresp(:,70+ttt);csize = size(all_resp_T,2)+1;
end
for ttt = 1:10
  all_resp_T(:,csize) = allresp(:,80+ttt);csize = size(all_resp_T,2)+1;
  all_resp_T(:,csize) = allresp(:,90+ttt);csize = size(all_resp_T,2)+1;
  all_resp_T(:,csize) = allresp(:,100+ttt);csize = size(all_resp_T,2)+1;
  all_resp_T(:,csize) = allresp(:,110+ttt);csize = size(all_resp_T,2)+1;
end
for ttt = 1:10
  all_resp_T(:,csize) = allresp(:,120+ttt);csize = size(all_resp_T,2)+1;
  all_resp_T(:,csize) = allresp(:,130+ttt);csize = size(all_resp_T,2)+1;
  all_resp_T(:,csize) = allresp(:,140+ttt);csize = size(all_resp_T,2)+1;
  all_resp_T(:,csize) = allresp(:,150+ttt);csize = size(all_resp_T,2)+1;
end
all_resp_T = [all_resp_T allresp(:,161:170)];csize = size(all_resp_T,2)+1;
for ttt = 1:10
  all_resp_T(:,csize) = allresp(:,170+ttt);csize = size(all_resp_T,2)+1;
  all_resp_T(:,csize) = allresp(:,180+ttt);csize = size(all_resp_T,2)+1;
  all_resp_T(:,csize) = allresp(:,190+ttt);csize = size(all_resp_T,2)+1;
end

%% Reviewer 1 last comment related
an = 1; pl = 1;

maskFrame = make_resp_frame(ei{an},pl,all_resp_T(an,:));
%%
%%
  ntrials = 50; %
    an = 4; pl = 1;

    sic = {[Lb Lbs];[Ab_On Abs_On];[Ab_Off Abs_Off];[Ab_Offc Abs_Offc];[Ar_On ArL_On Ars_On];[Ar_Off ArL_Off Ars_Off];[Ar_Offc ArL_Offc Ars_Offc]};
    sic = {[Lb];[Ab_On];[Ab_Off];[Ab_Offc];[Ar_On];[Ar_Off];[Ar_Offc]};
%     sic = {[Lb Lbs Ab_On Abs_On Ab_Off Abs_Off Ab_Offc Abs_Offc];[Ar_On ArL_On Ars_On Ar_Off ArL_Off Ars_Off Ar_Offc ArL_Offc Ars_Offc]};
%     sic = {[Ab_On Abs_On];[Ar_On ArL_On Ars_On]};
%     sic = {[Lb Lbs];[Ab_On Abs_On];[Ar_On ArL_On Ars_On];};
    clear all_gFR all_exc all_inh all_gV 
    prop_names = {'resp','N_Resp_Trials','zMI','zMINaN','HaFD','HiFD','cells_pooled'};
    cell_sel = {'good_FR','vals','exc','inh'};
    varName = {'all_gFR','all_gV','all_exc','all_inh'};
    for cii = 1:length(cell_sel);
        cmdTxt = sprintf('clear %s',varName{cii});eval(cmdTxt);
    end
    pni = 7;
    all_exc_inh = [];
    for ii = 1:length(sic)
        sit = sic{ii};
        tRs = o.Rs(:,sit);
        props1 = get_props_Rs(tRs,ntrials);
        for cii = 1:length(cell_sel)
            cmdTxt = sprintf('gFR = props1.%s;',cell_sel{cii});eval(cmdTxt);
            if pni == 1
                rf = find_percent(gFR);
                cmdTxt = sprintf('%s(:,ii) = mean(rf,2);',varName{cii}); eval(cmdTxt)
            else
                if pni == 7
                    cmdTxt = sprintf('%s(:,ii) = cell_list_op(gFR,[],''or'',1);',varName{cii});eval(cmdTxt);
                else
                    if pni == 3 || pni == 4
                        cmdtxt = sprintf('rf = props1.%s;',prop_names{3}); eval(cmdtxt);
                        if pni == 3
                            cmdTxt = sprintf('%s(:,ii) = mean(exec_fun_on_cell_mat(rf,''nanmean'',gFR),2);',varName{cii});eval(cmdTxt);
                        else
                            for rrr = 1:size(rf,1)
                                for ccc = 1:size(rf,2)
                                    temp = isnan(rf{rrr,ccc});
                                    rf1(rrr,ccc) = 100*sum(temp(gFR{rrr,ccc}))/size(rf{rrr,ccc},1);
                                end
                            end
                            cmdTxt = sprintf('%s(:,ii) = mean(rf1,2);',varName{cii}); eval(cmdTxt)
                        end
                    else
                        cmdtxt = sprintf('rf = props1.%s;',prop_names{pni}); eval(cmdtxt);
                        cmdTxt = sprintf('%s(:,ii) = mean(exec_fun_on_cell_mat(rf,''mean'',gFR),2);',varName{cii});eval(cmdTxt);
                    end
                end
            end
        end
        all_exc_inh = [all_exc_inh all_exc(:,ii) all_inh(:,ii)];
    end

% maskFrame = make_resp_frame(ei{an},pl,all_gFR(an,:));
% all_gFRs = all_gFR(an,1:2);
% maskFrame = make_resp_frame(ei{an},pl,all_gFRs);
% maskFrame = com_resp_frame(ei{an},pl,all_gFRs);
[maskFrame,all_masks,mimg] = make_resp_frame(ei{an},pl,all_inh(an,:));
% save_mean_img(ei);


maskFrame(maskFrame==0) = NaN;
% maskFrame = sum(all_masks,3);

ff = makeFigureRowsCols(107,[10 3 2.5 2.5],'RowsCols',[1 1],'spaceRowsCols',[0.2 0.2],'rightUpShifts',[0.15 0.13],'widthHeightAdjustment',[-10 -300]);
    set(gcf,'color','w');     ylims = [0 1];
    stp = 0.3; widths = [2 2 2 2 0.4 0.4]-0.2; gap = 1.75; adjust_axes(ff,ylims,stp,widths,gap,{'Euclidean Distance'});
    gap1 = 0.4; widths1 = 0.75;
mm = 0; 
MM = 7; %MM = max(maskFrame(:));
% rimg = imfuse(mimg,maskFrame);
h = imagesc(maskFrame);
axis equal
axis off
box off
% colorbar
[hc,hca] = putColorBar(gca,[0.09 0.07 -0.11 -0.15],{sprintf('%d',mm),sprintf('%d',MM)},6,'eastoutside',[0.07 0.05 0.07 0.05]);
  colormap jet
save_pdf(ff.hf,mData.pdf_folder,sprintf('topography.pdf'),600);
%%
ff = makeFigureRowsCols(107,[10 3 2.5 2.5],'RowsCols',[1 1],'spaceRowsCols',[0.2 0.2],'rightUpShifts',[0.15 0.13],'widthHeightAdjustment',[-10 -300]);
    set(gcf,'color','w');     ylims = [0 1];
    stp = 0.3; widths = [2 2 2 2 0.4 0.4]-0.2; gap = 1.75; adjust_axes(ff,ylims,stp,widths,gap,{'Euclidean Distance'});
    gap1 = 0.4; widths1 = 0.75;
mm = 0; 
MM = 7; %MM = max(maskFrame(:));
% rimg = imfuse(mimg,maskFrame);
h = imagesc(maskFrame);hold on;
axis equal
axis off
box off
nbins = 8; sz = size(maskFrame,1); bsize = size(maskFrame,1)/nbins;
sbin = 1:bsize:sz; ebin = bsize:bsize:sz;
for rr = 1:length(sbin)
    plot([1 512],[sbin(rr) sbin(rr)],'w');
    plot([sbin(rr) sbin(rr)],[1 512],'w');
end
% colorbar
[hc,hca] = putColorBar(gca,[0.09 0.07 -0.11 -0.15],{sprintf('%d',mm),sprintf('%d',MM)},6,'eastoutside',[0.07 0.05 0.07 0.05]);
  colormap jet
save_pdf(ff.hf,mData.pdf_folder,sprintf('topography.pdf'),600);
%%
tgv = all_gV(an,:);
for ii = 1:length(tgv)
    tgvm(:,ii) = tgv{1,ii};
end

%% quantification of salt and pepper distribution
% the idea is draw a circle around each cell and then see what is the id
% of the responsiveness from other sensorimotor events

all_nbins = [2 4 8 16 32];
for  an = 1:5 
      pl = 1;
        sic = {[Lb Lbs];[Ab_On Abs_On];[Ab_Off Abs_Off];[Ab_Offc Abs_Offc];[Ar_On ArL_On Ars_On];[Ar_Off ArL_Off Ars_Off];[Ar_Offc ArL_Offc Ars_Offc]};
        sic = {[Lb];[Ab_On];[Ab_Off];[Ab_Offc];[Ar_On];[Ar_Off];[Ar_Offc]};
    %     sic = {[Lb Lbs Ab_On Abs_On Ab_Off Abs_Off Ab_Offc Abs_Offc];[Ar_On ArL_On Ars_On Ar_Off ArL_Off Ars_Off Ar_Offc ArL_Offc Ars_Offc]};
    %     sic = {[Ab_On Abs_On];[Ar_On ArL_On Ars_On]};
    %     sic = {[Lb Lbs];[Ab_On Abs_On];[Ar_On ArL_On Ars_On];};
        clear all_gFR all_exc all_inh all_gV 
        prop_names = {'resp','N_Resp_Trials','zMI','zMINaN','HaFD','HiFD','cells_pooled'};
        cell_sel = {'good_FR','vals','exc','inh'};
        varName = {'all_gFR','all_gV','all_exc','all_inh'};
        for cii = 1:length(cell_sel);
            cmdTxt = sprintf('clear %s',varName{cii});eval(cmdTxt);
        end
        pni = 7;
        all_exc_inh = [];
        for ii = 1:length(sic)
            sit = sic{ii};
            tRs = o.Rs(:,sit);
            props1 = get_props_Rs(tRs,ntrials);
            for cii = 1:length(cell_sel)
                cmdTxt = sprintf('gFR = props1.%s;',cell_sel{cii});eval(cmdTxt);
                if pni == 1
                    rf = find_percent(gFR);
                    cmdTxt = sprintf('%s(:,ii) = mean(rf,2);',varName{cii}); eval(cmdTxt)
                else
                    if pni == 7
                        cmdTxt = sprintf('%s(:,ii) = cell_list_op(gFR,[],''or'',1);',varName{cii});eval(cmdTxt);
                    else
                        if pni == 3 || pni == 4
                            cmdtxt = sprintf('rf = props1.%s;',prop_names{3}); eval(cmdtxt);
                            if pni == 3
                                cmdTxt = sprintf('%s(:,ii) = mean(exec_fun_on_cell_mat(rf,''nanmean'',gFR),2);',varName{cii});eval(cmdTxt);
                            else
                                for rrr = 1:size(rf,1)
                                    for ccc = 1:size(rf,2)
                                        temp = isnan(rf{rrr,ccc});
                                        rf1(rrr,ccc) = 100*sum(temp(gFR{rrr,ccc}))/size(rf{rrr,ccc},1);
                                    end
                                end
                                cmdTxt = sprintf('%s(:,ii) = mean(rf1,2);',varName{cii}); eval(cmdTxt)
                            end
                        else
                            cmdtxt = sprintf('rf = props1.%s;',prop_names{pni}); eval(cmdtxt);
                            cmdTxt = sprintf('%s(:,ii) = mean(exec_fun_on_cell_mat(rf,''mean'',gFR),2);',varName{cii});eval(cmdTxt);
                        end
                    end
                end
            end
            all_exc_inh = [all_exc_inh all_exc(:,ii) all_inh(:,ii)];
        end

    % maskFrame = make_resp_frame(ei{an},pl,all_gFR(an,:));
    % all_gFRs = all_gFR(an,1:2);
    % maskFrame = make_resp_frame(ei{an},pl,all_gFRs);
    % maskFrame = com_resp_frame(ei{an},pl,all_gFRs);
    [maskFrame,all_masks,mimg] = make_resp_frame(ei{an},pl,all_inh(an,:));
    for bb = 1:length(all_nbins)
        nbins = all_nbins(bb); sz = size(maskFrame,1); bsize = size(maskFrame,1)/nbins;
        sbin = 1:bsize:sz; ebin = bsize:bsize:sz;
        all_uvals = []; uvals = [];
        for rr = 1:length(sbin)
            for cc = 1:length(sbin)
                uvals1 = 0;
                for mm = 1:size(all_masks,3)
                    subimage = all_masks(sbin(rr):ebin(rr),sbin(cc):ebin(cc),mm);
                    if sum(subimage(:))>0
                        uvals1 = uvals1 + 1;
                    end
                end
                uvals(rr,cc) = uvals1;
            end
        end
        all_uvals(an,:) = uvals(:);
        pra1(an,bb) = sum(all_uvals(an,:)>1)/size(all_uvals,2);
  end
%     for  an = 1:5
%     end
end
  %% make graph salt and pepper
  ff = makeFigureRowsCols(107,[10 3 2.5 2.5],'RowsCols',[1 1],'spaceRowsCols',[0.2 0.2],'rightUpShifts',[0.2 0.2],'widthHeightAdjustment',[-10 -300]);
    set(gcf,'color','w');     ylims = [0 1];
    stp = 0.55; widths = [2 2 2 2 0.4 0.4]-0.2; gap = 1.75; adjust_axes(ff,ylims,stp,widths,gap,{'Regions (%)'});
    gap1 = 0.4; widths1 = 0.75;
    [mVar,semVar] = findMeanAndStandardError(pra1);
    xd = all_nbins.^2;
    plot(xd,mVar,'o');
    errorbar(xd,mVar,semVar,'o');
    xlim([1 2000]);
    set(gca,'XScale','log','XTick',xd);
    xlabel('Number of Regions');ylabel({'Percentage of regions with more' 'than one type of response'});
    format_axes(gca);
    save_pdf(ff.hf,mData.pdf_folder,sprintf('topographyQ.pdf'),600);

%% Reviewer 1 last comment related

filename = fullfile(mData.pdf_folder,sprintf('movie'));
v = VideoWriter(filename,'MPEG-4');
open(v);
an = 4; pl = 1;
hf = figure(1000);clf;
ha = axes;
for ii = 1:200
  ha = showCells(ha,ei{an},pl,all_resp_T{an,ii},[0.3 0.3]);
  Fr = getframe(ha);
  writeVideo(v,Fr);
  writeVideo(v,Fr);
  writeVideo(v,Fr);
  pause(0.1);
end
close(v)


%% 1 off diagnoal (uniqe between adjacent trials) for responsive cells percentage graphs
while 1
    %%
    respV = [];
    for an = 1:5
        respV(an,:) = diag(all_CI_mat(:,:,an));
    end
    %% run one grand test
    [within,dvn,xlabels] = make_within_table({'cond','Trials'},[20,10]);
    dataT = make_between_table({respV},dvn);
    ra = RMA(dataT,within);
    ra.ranova;
    print_for_manuscript(ra)
    %% run one grand test for Reviewer 3
    red_respV = respV(:,[11:40 41:50 61:80]);
    red_respV = respV(:,[171:200 41:50 61:80]);
    [within,dvn,xlabels] = make_within_table({'BNB','Event','Trials'},[2,3,10]);
    dataT = make_between_table({red_respV},dvn);
    ra = RMA(dataT,within);
    ra.ranova;
    print_for_manuscript(ra)
    %% run 20 individual tests for each condition separately
    sts = 1:10:200; ses = 10:10:200;
    [within,dvn,xlabels] = make_within_table({'Trials'},[10]);
    for cni = 1:length(siG)
        t_respV = respV(:,sts(cni):ses(cni));
        dataT = make_between_table({t_respV},dvn);
        raR{cni} = RMA(dataT,within);
        raR{cni}.ranova;
        print_for_manuscript(raR{cni})
    end
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'cond','hsd'},[1.5 1 1]);
%     h(h==1) = 0;
	xdata = make_xdata([20],[1 1.5]);
    hf = get_figure(5,[8 7 3.5 1.25]);
    % s = generate_shades(length(bins)-1);
    tcolors = mData.dcolors;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',3,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.15);
%     set(hbs(1),'FaceColor',tcolors{1},'FaceAlpha',0.75); set(hbs(3),'FaceColor',tcolors{2},'FaceAlpha',0.75);
    maxY = maxY + 1;
    ylims = ylim;
    format_axes(gca);
%     htxt = text(0,maxY+6,sprintf('Any Condition (%d\x00B1%d%%),   All Conditions (%d%%)',round(mra),round(semra),round(mrall)),'FontSize',6);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[5 maxY]); format_axes(gca);
    xticks = xdata; 
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[5 15 25]); xtickangle(45)
    changePosition(gca,[0.0 0.01 0.05 0]); put_axes_labels(gca,{[],[0 0 0]},{{'Estimated','Marginal Means'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('responsivity.pdf'),600);
    
    %%
    mrespV = mean(respV); semrespV = std(respV)./sqrt(5);
    respTW = reshape(mrespV,10,20); mrespAct = respTW';
    respTW = reshape(semrespV,10,20); semrespAct = respTW';
    mrespActL = NaN;  semrespActL = NaN; 
    xtl = {''}; trialsStr = cellfun(@num2str,trials','un',0);
    xaxL = NaN; xticks = []; xtickL =[];
    for ii = 1:20
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
    plot(xlim,[nanmean(mrespActL) nanmean(mrespActL)],'m');
    iii=1;
    theinds = find(isnan(mrespActL));
    for ii = find(isnan(mrespActL))
        plot([ii ii],[11 26],'b-');
        if iii <= 20
        text(ii+2,29,sprintf('%s',xticklabels{iii}),'FontSize',6);
        indsS = (theinds(iii)+1):(theinds(iii+1)-1);
        shadedErrorBar(indsS,mrespActL(indsS),semrespActL(indsS),{},0.5)
        iii=iii+1;
        end
    end
    xlim([0 length(mrespActL)+1]); ylim([10 30]);
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

%% anova running
while 1
    %%
    aVar = [];
    for an = 1:5
        tvar = conjV(an,:);
        tvar(10:10:199) = NaN; tvar(isnan(tvar)) = []; 
        aVar(an,:) = tvar;
    end
    aVar_conj = aVar;
    %%
    [within,dvn,xlabels] = make_within_table({'cond','TrialsP'},[20,9]);
    dataT = make_between_table({aVar},dvn);
    rac = RMA(dataT,within);
    rac.ranova
     %% run one grand test for Reviewer 3
     red_aVar = aVar(:,[10:18 19:27 28:36 37:45 55:63 64:72]);
    [within,dvn,xlabels] = make_within_table({'BNB','Event','TrialsP'},[2,3,9]);
    dataT = make_between_table({red_aVar},dvn);
    racr = RMA(dataT,within);
    racr.ranova
    print_for_manuscript(racr)
     %% run 20 individual tests for each condition separately
    sts = 1:9:180; ses = 9:9:180;
    [within,dvn,xlabels] = make_within_table({'TrialsP'},[9]);
    for cni = 1:length(siG)
        t_aVar = aVar(:,sts(cni):ses(cni));
        dataT = make_between_table({t_aVar},dvn);
        raC{cni} = RMA(dataT,within);
        raC{cni}.ranova;
        print_for_manuscript(raC{cni})
    end
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,rac,{'cond','hsd'},[1.5 1 1]);
%     h(h==1) = 0;
	xdata = make_xdata([20],[1 1.5]);
    hf = get_figure(5,[8 7 3 1]);
    hf = get_figure(5,[8 7 3.5 1.25]);
    % s = generate_shades(length(bins)-1);
    tcolors = mData.dcolors;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',3,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.15);
%     set(hbs(1),'FaceColor',tcolors{1},'FaceAlpha',0.75); set(hbs(3),'FaceColor',tcolors{2},'FaceAlpha',0.75);
    maxY = maxY + 1;
    ylims = ylim;
    format_axes(gca);
%     htxt = text(0,maxY+6,sprintf('Any Condition (%d\x00B1%d%%),   All Conditions (%d%%)',round(mra),round(semra),round(mrall)),'FontSize',6);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[0 maxY]); format_axes(gca);
    xticks = xdata; 
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[10 20 30]); xtickangle(45)
    changePosition(gca,[0.01 0.01 0.05 0]); put_axes_labels(gca,{[],[0 0 0]},{{'Estimated','Marginal Means'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('conj_conds.pdf'),600);
    
    %%
    aVar = [];
    for an = 1:5
        tvar = comp1V(an,:);
        tvar(10:10:199) = NaN; tvar(isnan(tvar)) = []; 
        aVar(an,:) = tvar;
    end
    aVar_comp1 = aVar;
    %%
    [within,dvn,xlabels] = make_within_table({'cond','TrialsP'},[20,9]);
    dataT = make_between_table({aVar},dvn);
    rac1 = RMA(dataT,within);
    rac1.ranova
    %% run one grand test for Reviewer 3
     red_aVar = aVar(:,[10:18 19:27 28:36 37:45 55:63 64:72]);
    [within,dvn,xlabels] = make_within_table({'BNB','Event','TrialsP'},[2,3,9]);
    dataT = make_between_table({red_aVar},dvn);
    rac1r = RMA(dataT,within);
    rac1r.ranova
    print_for_manuscript(rac1r)
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,rac1r,{'BNB_by_Event','hsd'},[1.5 1 1]);
%     h(h==1) = 0;
	xdata = make_xdata([3 3],[1 1.5]);
    hf = get_figure(5,[8 7 3 1]);
    hf = get_figure(5,[8 7 2.5 1.25]);
    % s = generate_shades(length(bins)-1);
    tcolors = mData.dcolors;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',3,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.15);
%     set(hbs(1),'FaceColor',tcolors{1},'FaceAlpha',0.75); set(hbs(3),'FaceColor',tcolors{2},'FaceAlpha',0.75);
    maxY = maxY + 1;
    ylims = ylim;
    format_axes(gca);
%     htxt = text(0,maxY+6,sprintf('Any Condition (%d\x00B1%d%%),   All Conditions (%d%%)',round(mra),round(semra),round(mrall)),'FontSize',6);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[0 maxY]); format_axes(gca);
    xticks = xdata; 
    set(gca,'xtick',xticks,'xticklabels',{'AOn','AOff','Arb'},'ytick',[10 20 ]); xtickangle(45)
    set_bar_graph_sub_xtick_text(hf,gca,hbs,3,{'B','NB'});
    set_axes_top_text_no_line(hf,gca,'Complementary 1',[0.2 0.1 0 0]);
    changePosition(gca,[0.1 0.01 -0.07 0]); put_axes_labels(gca,{[],[0 0 0]},{{'Percentage','of Cells'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('comp1_conds.pdf'),600);
    
     %% run 20 individual tests for each condition separately
    sts = 1:9:180; ses = 9:9:180;
    [within,dvn,xlabels] = make_within_table({'TrialsP'},[9]);
    for cni = 1:length(siG)
        t_aVar = aVar(:,sts(cni):ses(cni));
        dataT = make_between_table({t_aVar},dvn);
        raC1{cni} = RMA(dataT,within);
        raC1{cni}.ranova;
        print_for_manuscript(raC1{cni})
    end
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,rac1,{'cond','hsd'},[1.5 1 1]);
%     h(h==1) = 0;
	xdata = make_xdata([20],[1 1.5]);
    hf = get_figure(5,[8 7 3 1]);
    hf = get_figure(5,[8 7 3.5 1.25]);
    % s = generate_shades(length(bins)-1);
    tcolors = mData.dcolors;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',3,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.15);
%     set(hbs(1),'FaceColor',tcolors{1},'FaceAlpha',0.75); set(hbs(3),'FaceColor',tcolors{2},'FaceAlpha',0.75);
    maxY = maxY + 1;
    ylims = ylim;
    format_axes(gca);
%     htxt = text(0,maxY+6,sprintf('Any Condition (%d\x00B1%d%%),   All Conditions (%d%%)',round(mra),round(semra),round(mrall)),'FontSize',6);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[0 maxY]); format_axes(gca);
    xticks = xdata; 
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[10 20 30]); xtickangle(45)
    changePosition(gca,[0.01 0.01 0.05 0]); put_axes_labels(gca,{[],[0 0 0]},{{'Estimated','Marginal Means'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('comp1_conds.pdf'),600);
    %%
    aVar = [];
    for an = 1:5
        tvar = comp2V(an,:);
        tvar(10:10:199) = NaN; tvar(isnan(tvar)) = []; 
        aVar(an,:) = tvar;
    end
    aVar_comp2 = aVar;
    %%
    [within,dvn,xlabels] = make_within_table({'cond','TrialsP'},[20,9]);
    dataT = make_between_table({aVar},dvn);
    rac2 = RMA(dataT,within);
    rac2.ranova
    %% run one grand test for Reviewer 3
     red_aVar = aVar(:,[10:18 19:27 28:36 37:45 55:63 64:72]);
    [within,dvn,xlabels] = make_within_table({'BNB','Event','TrialsP'},[2,3,9]);
    dataT = make_between_table({red_aVar},dvn);
    rac2r = RMA(dataT,within);
    rac2r.ranova
    print_for_manuscript(rac2r)
     %% run 20 individual tests for each condition separately
    sts = 1:9:180; ses = 9:9:180;
    [within,dvn,xlabels] = make_within_table({'TrialsP'},[9]);
    for cni = 1:length(siG)
        t_aVar = aVar(:,sts(cni):ses(cni));
        dataT = make_between_table({t_aVar},dvn);
        raC2{cni} = RMA(dataT,within);
        raC2{cni}.ranova;
        print_for_manuscript(raC2{cni})
    end
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,rac2r,{'BNB','hsd'},[1.5 1 1]);
%     h(h==1) = 0;
	xdata = make_xdata([2],[1 1.5]);
    hf = get_figure(5,[8 7 3 1]);
    hf = get_figure(5,[8 7 1.5 1.25]);
    % s = generate_shades(length(bins)-1);
    tcolors = mData.dcolors;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',3,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.15);
%     set(hbs(1),'FaceColor',tcolors{1},'FaceAlpha',0.75); set(hbs(3),'FaceColor',tcolors{2},'FaceAlpha',0.75);
    maxY = maxY + 1;
    ylims = ylim;
    format_axes(gca);
%     htxt = text(0,maxY+6,sprintf('Any Condition (%d\x00B1%d%%),   All Conditions (%d%%)',round(mra),round(semra),round(mrall)),'FontSize',6);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[0 maxY]); format_axes(gca);
    xticks = xdata; 
    set(gca,'xtick',xticks,'xticklabels',{'B','NB'},'ytick',[10 20 ]); xtickangle(45)
    changePosition(gca,[0.15 0.01 -0.07 0]); put_axes_labels(gca,{[],[0 0 0]},{{'Percentage','of Cells'},[0 0 0]});
    set_axes_top_text_no_line(hf,gca,'Complementary 2',[0.2 0.1 0 0]);
    save_pdf(hf,mData.pdf_folder,sprintf('comp2_conds.pdf'),600);
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,rac2,{'TrialsP','hsd'},[1.5 1 1]);
%     h(h==1) = 0;
	xdata = make_xdata([9],[1 1.5]);
    hf = get_figure(5,[8 7 2 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = mData.colors;
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
    xtltp = {'1-2','2-3','3-4','4-5','5-6','6-7','7-8','8-9','9-10'};
    set(gca,'xtick',xticks,'xticklabels',xtltp,'ytick',[5 10 20 30]); xtickangle(45)
    changePosition(gca,[0.02 0.11 0.05 -0.1]); put_axes_labels(gca,{'Trial-Pairs',[0 0 0]},{{'Cells (%)'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('comp2_trials.pdf'),600);
    %%
    break;
end


%% plot reviewer 3 point 3
eN = [2 5];
sss = 1:9:180; eee = 9:9:180;
aVar = [aVar_conj(:,[sss(eN(1)):eee(eN(1)) sss(eN(2)):eee(eN(2))]) aVar_comp1(:,[sss(eN(1)):eee(eN(1)) sss(eN(2)):eee(eN(2))]) aVar_comp2(:,[sss(eN(1)):eee(eN(1)) sss(eN(2)):eee(eN(2))])];

[within,dvn,xlabels] = make_within_table({'BNB','CT','TrialsP'},[2,3,9]);
dataT = make_between_table({aVar},dvn);
raccc = RMA(dataT,within);
raccc.ranova
print_for_manuscript(raccc)
%%
magfac = mData.magfac;
ff = makeFigureRowsCols(108,[5 5 1.25 1],'RowsCols',[1 1],'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.01 0.25],'widthHeightAdjustment',[10 -320]);
ff = makeFigureRowsCols(108,[5 5 2.25 1],'RowsCols',[1 1],'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.01 0.25],'widthHeightAdjustment',[10 -320]);
MY = 20; ysp = 3; mY = 0; titletxt = 'Responsivity'; ylabeltxt = {'Percent of Cells'};

stp = 0.3;magfac; widths = ([0.75 ])*magfac; gap = 0.01*magfac;
stp = 0.3;magfac; widths = ([1.75 ])*magfac; gap = 0.01*magfac;
adjust_axes(ff,[mY MY],stp,widths,gap,{''});

[xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,raccc,{'TrialsP','hsd'},[1.5 1 1]);
xdata = make_xdata([2],[1 1.5]);   
xdata = make_xdata([3],[1 1.5]);   
xdata = make_xdata([9],[1 1.5]); 
%     mVar = ra.est_marginal_means.Mean; semVar = ra.est_marginal_means.Formula_StdErr;
% tcolors = [temp_tcolors temp_tcolors];
%     combs = ra.mcs.combs; p = ra.mcs.p; h = p<0.00005;
%     xdata = [1:(length(mVar)/2) ((length(mVar)/2)+1+(1:(length(mVar)/2)))];
tcolors = mData.dcolors(4:end);
%     hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 6.9 2],'color','w');
axes(ff.h_axes(1,1))
[hbs,maxY]  = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
    'ySpacing',4,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
    'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.5,'sigLinesStartYFactor',0.05);
% make_bars_hollow(hbs(13:end));
set_axes_limits(gca,[0.25 xdata(end)+.75],[mY MY]); format_axes(gca); xticks = xdata; 
set(gca,'xtick',xticks,'xticklabels',{'B','NB'});xtickangle(45)
set(gca,'xtick',xticks,'xticklabels',{'Conj','Comp1','Comp2'});xtickangle(35)
set(gca,'xtick',xticks,'xticklabels',{'T12','T23','T34','T45','T56','T67','T78','T89','T9T'});xtickangle(45)
put_axes_labels(gca,{[],[0 0 0]},{ylabeltxt,[0 0 0]});
ht = set_axes_top_text_no_line(gcf,gca,'',[0 -0.051 0 0]);
save_pdf(ff.hf,mData.pdf_folder,sprintf('conjcompcomp.pdf'),600);
%     changePosition(gca,[-0.03 0.03 0.11 -0.01]);
% set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,3,{'C1','C2','C3','C4','C1','C2','C3','C4'},{[-0.01 0.02]});
% set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,12,{'Control','APP'},{[-0.12 0]});

%%
eN = [3 7];
sss = 1:9:180; eee = 9:9:180;
aVar = [aVar_conj(:,[sss(eN(1)):eee(eN(1)) sss(eN(2)):eee(eN(2))]) aVar_comp1(:,[sss(eN(1)):eee(eN(1)) sss(eN(2)):eee(eN(2))]) aVar_comp2(:,[sss(eN(1)):eee(eN(1)) sss(eN(2)):eee(eN(2))])];

[within,dvn,xlabels] = make_within_table({'BNB','CT','TrialsP'},[2,3,9]);
dataT = make_between_table({aVar},dvn);
raccc = RMA(dataT,within);
raccc.ranova
print_for_manuscript(raccc)
%%
magfac = mData.magfac;
ff = makeFigureRowsCols(108,[5 5 2.25 1],'RowsCols',[1 1],'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.01 0.35],'widthHeightAdjustment',[10 -420]);
MY = 35; ysp = 3; mY = 0; titletxt = 'Responsivity'; ylabeltxt = {'Percent of Cells'};

stp = 0.3;magfac; widths = ([1.75])*magfac; gap = 0.01*magfac;
adjust_axes(ff,[mY MY],stp,widths,gap,{''});

[xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,raccc,{'BNB_by_CT','hsd'},[1.5 1 1]);
xdata = make_xdata([3 3],[1 1.5]);   

%     mVar = ra.est_marginal_means.Mean; semVar = ra.est_marginal_means.Formula_StdErr;
% tcolors = [temp_tcolors temp_tcolors];
%     combs = ra.mcs.combs; p = ra.mcs.p; h = p<0.00005;
%     xdata = [1:(length(mVar)/2) ((length(mVar)/2)+1+(1:(length(mVar)/2)))];
tcolors = [mData.dcolors(1:3) mData.dcolors(1:3)];
%     hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 6.9 2],'color','w');
axes(ff.h_axes(1,1))
[hbs,maxY]  = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
    'ySpacing',4,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
    'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.5,'sigLinesStartYFactor',0.05);
% make_bars_hollow(hbs(13:end));
set_axes_limits(gca,[0.25 xdata(end)+.75],[mY MY]); format_axes(gca); xticks = xdata; 
set(gca,'xtick',xticks,'xticklabels',{'B','NB'});xtickangle(45)
set(gca,'xtick',xticks,'xticklabels',{'Conj','Comp1','Comp2'});xtickangle(30)
put_axes_labels(gca,{[],[0 0 0]},{ylabeltxt,[0 0 0]});
ht = set_axes_top_text_no_line(gcf,gca,'',[0 -0.051 0 0]);
set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,3,{'Brake','No-Brake'},{[0.01 0.02]});
save_pdf(ff.hf,mData.pdf_folder,sprintf('conjcompcomp.pdf'),600);
%%
eN = [4 8];
sss = 1:9:180; eee = 9:9:180;
aVar = [aVar_conj(:,[sss(eN(1)):eee(eN(1)) sss(eN(2)):eee(eN(2))]) aVar_comp1(:,[sss(eN(1)):eee(eN(1)) sss(eN(2)):eee(eN(2))]) aVar_comp2(:,[sss(eN(1)):eee(eN(1)) sss(eN(2)):eee(eN(2))])];

[within,dvn,xlabels] = make_within_table({'BNB','CT','TrialsP'},[2,3,9]);
dataT = make_between_table({aVar},dvn);
raccc = RMA(dataT,within);
raccc.ranova
print_for_manuscript(raccc)
%%
magfac = mData.magfac;
ff = makeFigureRowsCols(108,[5 5 2.25 1],'RowsCols',[1 1],'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.01 0.35],'widthHeightAdjustment',[10 -420]);
MY = 30; ysp = 3; mY = 0; titletxt = 'Responsivity'; ylabeltxt = {'Percent of Cells'};

stp = 0.3;magfac; widths = ([1.75])*magfac; gap = 0.01*magfac;
adjust_axes(ff,[mY MY],stp,widths,gap,{''});

[xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,raccc,{'BNB_by_CT','hsd'},[1.5 1 1]);
xdata = make_xdata([3 3],[1 1.5]);   

%     mVar = ra.est_marginal_means.Mean; semVar = ra.est_marginal_means.Formula_StdErr;
% tcolors = [temp_tcolors temp_tcolors];
%     combs = ra.mcs.combs; p = ra.mcs.p; h = p<0.00005;
%     xdata = [1:(length(mVar)/2) ((length(mVar)/2)+1+(1:(length(mVar)/2)))];
tcolors = [mData.dcolors(1:3) mData.dcolors(1:3)];
%     hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 6.9 2],'color','w');
axes(ff.h_axes(1,1))
[hbs,maxY]  = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
    'ySpacing',4,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
    'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.5,'sigLinesStartYFactor',0.05);
% make_bars_hollow(hbs(13:end));
set_axes_limits(gca,[0.25 xdata(end)+.75],[mY MY]); format_axes(gca); xticks = xdata; 
set(gca,'xtick',xticks,'xticklabels',{'B','NB'});xtickangle(45)
set(gca,'xtick',xticks,'xticklabels',{'Conj','Comp1','Comp2'});xtickangle(30)
put_axes_labels(gca,{[],[0 0 0]},{ylabeltxt,[0 0 0]});
ht = set_axes_top_text_no_line(gcf,gca,'',[0 -0.051 0 0]);
set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,3,{'Brake','No-Brake'},{[0.01 0.02]});
save_pdf(ff.hf,mData.pdf_folder,sprintf('conjcompcomp.pdf'),600);