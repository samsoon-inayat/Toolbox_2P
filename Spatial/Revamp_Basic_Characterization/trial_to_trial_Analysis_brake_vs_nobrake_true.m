function trial_to_trial_Analysis_brake_vs_nobrake

%% find spatial trial to trial correlation
while 1
    trialNums = [1:10];
   si = [Lb Ab_On Ab_Off Ar_On Ar_Off ArL_On ArL_Off Ars_On Ars_Off Lbs Abs_On Abs_Off ArL_L];
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
    respG = avgProps.good_FR_and_untuned;
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

    [OIo,mOI,semOI,OI_mato,p_vals,h_vals,all_CI,mCI,semCI,all_CI_mat,uni] = get_overlap_index(allresp,0.5,0.05);
    mOI = mCI; semOI = semCI;
    disp('Done');
    %%
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
%     ff = makeFigureRowsCols(2020,[10 4 6 1.5],'RowsCols',[1 2],'spaceRowsCols',[0.1 0.01],'rightUpShifts',[0.05 0.13],'widthHeightAdjustment',[-240 -150]);
    hf = get_figure(6,[8 3 2.75 2.75]);
%     hf = get_figure(6,[8 3 4 4]);
    %
%     axes(ff.h_axes(1));
    im1 = imagesc(mOI,[minI,maxI]);    im1.AlphaData = imAlpha;
    
    for ii = 1:12
        plot([(10.5+((ii-1)*10)) (10.5+((ii-1)*10))],[0 130.5],'w','linewidth',0.1); 
        plot([0 130.5],[(10.5+((ii-1)*10)) (10.5+((ii-1)*10))],'w','linewidth',0.1); 
    end
%     for ii = [2 6 9 12]   
%         plot([(10.5+((ii-1)*10)) (10.5+((ii-1)*10))],[0 130.5],'k','linewidth',1); 
%         plot([0 130.5],[(10.5+((ii-1)*10)) (10.5+((ii-1)*10))],'k','linewidth',1); 
%     end
    set(gca,'color',0.5*[1 1 1]);    colormap parula;    %axis equal
    format_axes(gca);
    set_axes_limits(gca,[0.5 sz+0.5],[0.5 sz+0.5]);
    ttxl = rasterNamesTxt(siG);
    xtickvals = 5:10:130;%[5 15 25 60 100 115 125];
    xticklabels = rasterNamesTxt(siG);

    set(gca,'xtick',xtickvals,'ytick',xtickvals,'xticklabels',xticklabels,'yticklabels',xticklabels,'Ydir','normal'); xtickangle(45);%ytickangle(45);
    yyaxis right
    set(gca,'ytick',xtickvals,'yticklabels',yticklabels,'tickdir','out');
    box off
    changePosition(gca,[-0.01 0.00 0.037 0.033]);
    hc = putColorBar(gca,[0.1 -0.08 -0.2 0.03],{sprintf('%.1f',minI),sprintf('%.1f',maxI)},6,'northoutside',[0.07 0.09 0.02 0.09]);
    colormap jet
    save_pdf(hf,mData.pdf_folder,sprintf('OI_Map_mean_tu_spatial.pdf'),600);
    %%
    break;
end
%% along diagnol (responsiveness)
while 1
mask = diag(mCI);
hf = figure(200);clf;
set(hf,'units','inches','position',[5 5 3.5 1.5]);
plot(1:130,mask,'m');hold on;
mconj = nanmean(mask);
plot([1 130],[mconj mconj],'k');
xlim([1 130]);
for ii = 1:13
    plot([ii*10 ii*10],ylim,'b-');
%     text(ii*10-5,23.5,sprintf('%s',xticklabels{ii}),'FontSize',6);
end
xticks = 5:10:130;
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
respTW = reshape(mask,10,13);
respTW = respTW';
hf = figure(100);clf;
set(hf,'units','inches','position',[5 5 3 1.5]);
tcolors = mData.dcolors(1:13);
for ii = 1:13
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
mask(10:10:129) = NaN;
mask(130) = NaN;
mconj = nanmean(mask);
hf = figure(100);clf;
set(hf,'units','inches','position',[5 5 3.5 0.75]);
plot(1:130,mask,'m');hold on;
plot([1 130],[mconj mconj],'k');
xlim([1 130]);
for ii = 1:13
    plot([ii*10 ii*10],[4 11],'b--');
%     text(ii*10-5,23.5,sprintf('%s',xticklabels{ii}),'FontSize',6);
end
xticks = 5:10:130;
set(gca,'Xtick',xticks,'XTickLabels',xticklabels);xtickangle(30);
ylabel('Cells (%)');box off;
format_axes(gca);
save_pdf(hf,mData.pdf_folder,sprintf('tria_to_trial_conj.pdf'),600);
break;
end

%% 1 off diagnoal (conjunctive between adjacent trials) trial-pair wise
while 1
mask = diag(mCI,1);
mask(10:10:129) = NaN;
mask(130) = NaN;
mask(isnan(mask)) = [];
respTW = reshape(mask,9,13);
respTW = respTW';
hf = figure(100);clf;
set(hf,'units','inches','position',[5 5 3 1.5]);
tcolors = mData.dcolors(1:13);
for ii = 1:13
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
        inds = trialN:10:129;
        tvals(an,:) = vals(inds);
    end
    
    [within,dvn,xlabels] = make_within_table({'cond'},[13]);
    dataT = make_between_table({tvals},dvn);
    ra = RMA(dataT,within);
    ra.ranova;
    pvals(trialN) = ra.ranova.pValue_sel(3);
    eta2(trialN) = ra.eta2;
        
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'cond','hsd'},[1.5 1 1]);
%     h(h==1) = 0;
	xdata = make_xdata([13],[1 1.5]);
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
%     ff = makeFigureRowsCols(2020,[10 4 6 1.5],'RowsCols',[1 2],'spaceRowsCols',[0.1 0.01],'rightUpShifts',[0.05 0.13],'widthHeightAdjustment',[-240 -150]);
    hf = get_figure(6,[8 3 2.75 2.75]);
%     hf = get_figure(6,[8 3 4 4]);
    %
%     axes(ff.h_axes(1));
    im1 = imagesc(mOI,[minI,maxI]);    im1.AlphaData = imAlpha;
    
    for ii = 1:12
        plot([(10.5+((ii-1)*10)) (10.5+((ii-1)*10))],[0 130.5],'w','linewidth',0.1); 
        plot([0 130.5],[(10.5+((ii-1)*10)) (10.5+((ii-1)*10))],'w','linewidth',0.1); 
    end
%     for ii = [2 6 9 12]   
%         plot([(10.5+((ii-1)*10)) (10.5+((ii-1)*10))],[0 130.5],'k','linewidth',1); 
%         plot([0 130.5],[(10.5+((ii-1)*10)) (10.5+((ii-1)*10))],'k','linewidth',1); 
%     end
    set(gca,'color',0.5*[1 1 1]);    colormap parula;    %axis equal
    format_axes(gca);
    set_axes_limits(gca,[0.5 sz+0.5],[0.5 sz+0.5]);
    ttxl = rasterNamesTxt(siG);
    xtickvals = 5:10:130;%[5 15 25 60 100 115 125];
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
mask(10:10:129) = NaN;
mask(130) = NaN;
mconj = nanmean(mask);
hf = figure(100);clf;
set(hf,'units','inches','position',[5 5 3.5 0.75]);
plot(1:130,mask,'m');hold on;
plot([1 130],[mconj mconj],'k');
xlim([1 130]);
for ii = 1:13
    plot([ii*10 ii*10],ylim,'b--');
%     text(ii*10-5,23.5,sprintf('%s',xticklabels{ii}),'FontSize',6);
end
xticks = 5:10:130;
set(gca,'Xtick',xticks,'XTickLabels',xticklabels);xtickangle(30);
ylabel('Cells (%)');box off;
format_axes(gca);
save_pdf(hf,mData.pdf_folder,sprintf('tria_to_trial_unique.pdf'),600);
break;
end

%% 1 off diagnoal (unique ) trial-pair wise
while 1

mask = diag(mOI,-1);
mask(10:10:129) = NaN;
mask(130) = NaN;
mask(isnan(mask)) = [];
respTW = reshape(mask,9,13);
respTW = respTW';
hf = figure(100);clf;
set(hf,'units','inches','position',[5 5 3 1.5]);
tcolors = mData.dcolors(1:13);
for ii = 1:13
    plot(1:9,respTW(ii,:),'color',mData.dcolors{ii});hold on;
end
mconj = nanmean(respTW); semconj = std(respTW)./sqrt(5);
h = shadedErrorBar(1:9,mconj,semconj,{'color','k','linewidth',1.5},0.5);
xlim([1 11]); ylim([4,18]);
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
        inds = trialN:10:129;
        tvals(an,:) = vals(inds);
    end
    
    [within,dvn,xlabels] = make_within_table({'cond'},[13]);
    dataT = make_between_table({tvals},dvn);
    ra = RMA(dataT,within);
    ra.ranova;
    pvals(trialN) = ra.ranova.pValue_sel(3);
    eta2(trialN) = ra.eta2;
        
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'cond','hsd'},[1.5 1 1]);
%     h(h==1) = 0;
	xdata = make_xdata([13],[1 1.5]);
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
mask(10:10:129) = NaN;
mask(130) = NaN;
mconj = nanmean(mask);

masku = diag(mOI,-1);
masku(10:10:129) = NaN;
masku(130) = NaN;
mconju = nanmean(masku);

hf = figure(100);clf;
set(hf,'units','inches','position',[5 5 3.5 0.75]);
plot(1:130,mask,'m');hold on;
plot(1:130,masku,'c');
plot([1 130],[mconj mconj],'k');
xlim([1 130]);
for ii = 1:13
    plot([ii*10 ii*10],ylim,'b--');
%     text(ii*10-5,23.5,sprintf('%s',xticklabels{ii}),'FontSize',6);
end
xticks = 5:10:130;
set(gca,'Xtick',xticks,'XTickLabels',xticklabels);xtickangle(30);
ylabel('Cells (%)');box off;
format_axes(gca);
save_pdf(hf,mData.pdf_folder,sprintf('tria_to_trial_disrupted.pdf'),600);
break;
end

%% 1 off diagnoal (disrupted ) trial-pair wise
while 1

mask = diag(mOI,1);
mask(10:10:129) = NaN;
mask(130) = NaN;
mask(isnan(mask)) = [];
respTW = reshape(mask,9,13);
respTW = respTW';
hf = figure(100);clf;
set(hf,'units','inches','position',[5 5 3 1.5]);
tcolors = mData.dcolors(1:13);
for ii = 1:13
    plot(1:9,respTW(ii,:),'color',mData.dcolors{ii});hold on;
end
mconj = nanmean(respTW); semconj = std(respTW)./sqrt(5);
h = shadedErrorBar(1:9,mconj,semconj,{'color','k','linewidth',1.5},0.5);
xlim([1 11]); ylim([4,18]);
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
        inds = trialN:10:129;
        tvals(an,:) = vals(inds);
    end
    
    [within,dvn,xlabels] = make_within_table({'cond'},[13]);
    dataT = make_between_table({tvals},dvn);
    ra = RMA(dataT,within);
    ra.ranova;
    pvals(trialN) = ra.ranova.pValue_sel(3);
    eta2(trialN) = ra.eta2;
        
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'cond','hsd'},[1.5 1 1]);
%     h(h==1) = 0;
	xdata = make_xdata([13],[1 1.5]);
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
while 1
    respAct = diag(mCI);
mask = diag(mOI,-1);mask(10:10:129) = NaN;mask(130) = NaN;
maskd = diag(mOI,1);maskd(10:10:129) = NaN;maskd(130) = NaN;
maskc = diag(mCI,1);maskc(10:10:129) = NaN;maskc(130) = NaN;
mconj = nanmean(mask); mconjd = nanmean(maskd); mconjc = nanmean(maskc);

hf = figure(100);clf;
set(hf,'units','inches','position',[5 5 6.9 1.5]);
plot(1:130,mask,'m');hold on;
plot(1:130,respAct,'r');
plot(1:130,maskd,'c');
plot(1:130,maskc,'k');
plot([1 130],[mconj mconj],'m'); plot([1 130],[mconjd mconjd],'c'); plot([1 130],[mconjc mconjc],'k');
xlim([0 131]); ylim([2 25]);
for ii = 1:13
    plot([ii*10 ii*10],ylim,'b-');
%     text(ii*10-5,23.5,sprintf('%s',xticklabels{ii}),'FontSize',6);
end
xticks = sort([5:10:130 0:10:130]);
% set(gca,'Xtick',xticks,'XTickLabels',xticklabels);xtickangle(30);
set(gca,'Xtick',xticks)
ylabel('Cells (%)');box off;
format_axes(gca);
changePosition(gca,[-0.08 0.1 0.17 -0.1]);
save_pdf(hf,mData.pdf_folder,sprintf('tria_to_trial_unique.pdf'),600);
break;
end
