function trial_to_trial_Analysis_brake_vs_nobrake
%%
G = 'C'; o = oC;
G = 'A'; o = oA;
n = 0;
%% find spatial trial to trial correlation
while 1
    trialNums = [1:10];
   si = [C1_t_D C1_i_T C2_t_D C2_i_T C3_t_D C3_i_T C4_t_D C4_i_T];
   ind = 1;
   for rr = 1:4
       for cc = 1:10
           event_type{ind} = sprintf('D%d%d',rr,cc-1);
           ind = ind + 1;
       end
       for cc = 1:10
           event_type{ind} = sprintf('T%d%d',rr,cc-1);
           ind = ind + 1;
       end
   end
   
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
    
    [OIo,mOI,semOI,OI_mato,p_vals,h_vals,all_CI,mCI,semCI,all_CI_mat,uni] = get_overlap_index(allresp,0.5,0.05);
    
    mUni = nanmean(uni,3); mUni1 = tril(mUni,-1) + tril(mUni,-1)'; mUni2 = triu(mUni,1) + triu(mUni,1)'; mmUni = min(mUni(:)); MmUni = max(mUni(:));
    mOI = mCI; semOI = semCI;
    mSel = mCI;
    mOI = mSel;
    sz = size(mOI,1);
%     mOI = OI_mat(:,:,4);
    oM = ones(size(mOI));
    mask1 = (triu(oM,0) & tril(oM,0)); mOI(mask1==1) = NaN;
    maxI = 35;max([mOI(:);semOI(:)]);
    minI = min([mOI(:);semOI(:)]);
%     minI = 0; maxI = 0.6;
    mask = tril(NaN(size(mOI)),0); mask(mask==0) = 1; 
%     mOI = mOI .* mask;
    imAlpha=ones(size(mOI));    %imAlpha(isnan(mask))=0.25; 
    imAlpha(mask1 == 1) = 0;
%     ff = makeFigureRowsCols(2020,[10 4 6 1.5],'RowsCols',[1 2],'spaceRowsCols',[0.1 0.01],'rightUpShifts',[0.05 0.8],'widthHeightAdjustment',[-240 -150]);
    hf = get_figure(6,[8 3 3.5 3.5]);
%     hf = get_figure(6,[8 3 4 4]);
    %
%     axes(ff.h_axes(1));
    im1 = imagesc(mOI,[minI,maxI]);    im1.AlphaData = imAlpha;
    
    for ii = 1:18
        plot([(10.5+((ii-1)*10)) (10.5+((ii-1)*10))],[0 180.5],'w','linewidth',0.1); 
        plot([0 180.5],[(10.5+((ii-1)*10)) (10.5+((ii-1)*10))],'w','linewidth',0.1); 
    end
%     for ii = [2 6 9 12]   
%         plot([(10.5+((ii-1)*10)) (10.5+((ii-1)*10))],[0 80.5],'k','linewidth',1); 
%         plot([0 80.5],[(10.5+((ii-1)*10)) (10.5+((ii-1)*10))],'k','linewidth',1); 
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
    save_pdf(hf,mData.pdf_folder,sprintf('OI_Map_mean_tu_spatial_%c.pdf',G),600);
    
    
     %%
     mCI1 = mSel;
     mask1 = (triu(oM,0) & tril(oM,0)); mCI1(mask1==1) = NaN;
     sts = 1:10:size(mOI,2); ses = 10:10:size(mOI,2);
     clear mOI;
     for rr = 1:length(siG)
         for cc = 1:length(siG)
             mOI(rr,cc) = mean(mCI1(sts(rr):ses(rr),sts(cc):ses(cc)),'All');
         end
     end
%      mOI = mCI; semOI = semCI;
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
%     ff = makeFigureRowsCols(2020,[10 4 6 1.5],'RowsCols',[1 2],'spaceRowsCols',[0.1 0.01],'rightUpShifts',[0.05 0.8],'widthHeightAdjustment',[-240 -150]);
    hf = get_figure(6,[8 3 3.5 3.5]);
%     hf = get_figure(6,[8 3 4 4]);
    %
%     axes(ff.h_axes(1));
    im1 = imagesc(mOI,[minI,maxI]);    im1.AlphaData = imAlpha;
    
    for ii = 1:15
        plot([(10.5+((ii-1)*10)) (10.5+((ii-1)*10))],[0 80.5],'w','linewidth',0.1); 
        plot([0 80.5],[(10.5+((ii-1)*10)) (10.5+((ii-1)*10))],'w','linewidth',0.1); 
    end
%     for ii = [2 6 9 12]   
%         plot([(10.5+((ii-1)*10)) (10.5+((ii-1)*10))],[0 80.5],'k','linewidth',1); 
%         plot([0 80.5],[(10.5+((ii-1)*10)) (10.5+((ii-1)*10))],'k','linewidth',1); 
%     end
    set(gca,'color',0.5*[1 1 1]);    colormap parula;    %axis equal
    format_axes(gca);
    set_axes_limits(gca,[0.5 sz+0.5],[0.5 sz+0.5]);
    ttxl = rasterNamesTxt(siG);
    xtickvals = 5:10:150;%[5 15 25 60 100 115 125];
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
    mOI1 = mOI;
%     mOI1(isnan(mOI1)) = 1;
    Di = pdist(mOI1,@naneucdist);
%     tree = linkage(mOI1,'average','euclidean');
    tree = linkage(Di,'average');
    figure(hf);clf
    [H,T,TC] = dendrogram(tree,81,'Orientation','top','ColorThreshold','default');
    hf = gcf;
    txl = event_type;
    set(hf,'Position',[7 3 6.9 1.5]);
    set(H,'linewidth',1);
    set(gca,'xticklabels',txl(TC));xtickangle(45);
    format_axes(gca);
    hx = ylabel('Eucledian Distance');changePosition(hx,[0 -0.1 0]);
    changePosition(gca,[-0.05 0.0 0.09 0.05]);
    save_pdf(hf,mData.pdf_folder,sprintf('OI_Map_cluster_%c.pdf',G),600);
    %%
    break;
end


%% 1 off diagnoal (uniqe between adjacent trials)
while 1
    
    for an = 1:5
        respRV(an,:) = diag(all_CI_mat(:,:,an));
        conjV(an,:) = diag(all_CI_mat(:,:,an),1);
        comp1V(an,:) = diag(uni(:,:,an),1);
        comp2V(an,:) = diag(uni(:,:,an),-1);
    end
    respV = respRV;
    mrespV = mean(respV); semrespV = std(respV)./sqrt(5);
    respTW = reshape(mrespV,10,8); mrespAct = respTW'; 
    respTW = reshape(semrespV,10,8); semrespAct = respTW';
    
    respV = conjV;
    mrespV = mean(respV); semrespV = std(respV)./sqrt(5);
    mrespV(10:10:79) = NaN; mrespV(80) = NaN; mrespV(isnan(mrespV)) = []; 
    semrespV(10:10:79) = NaN; semrespV(80) = NaN; semrespV(isnan(semrespV)) = []; 
    respTW = reshape(mrespV,9,8); mconjAct = respTW'; 
    respTW = reshape(semrespV,9,8); semconjAct = respTW';
    
    respV = comp1V;
    mrespV = mean(respV); semrespV = std(respV)./sqrt(5);
    mrespV(10:10:79) = NaN; mrespV(80) = NaN; mrespV(isnan(mrespV)) = []; 
    semrespV(10:10:79) = NaN; semrespV(80) = NaN; semrespV(isnan(semrespV)) = []; 
    respTW = reshape(mrespV,9,8); mcomp1Act = respTW'; 
    respTW = reshape(semrespV,9,8); semcomp1Act = respTW';
    
    respV = comp2V;
    mrespV = mean(respV); semrespV = std(respV)./sqrt(5);
    mrespV(10:10:79) = NaN; mrespV(80) = NaN; mrespV(isnan(mrespV)) = []; 
    semrespV(10:10:79) = NaN; semrespV(80) = NaN; semrespV(isnan(semrespV)) = []; 
    respTW = reshape(mrespV,9,8); mcomp2Act = respTW'; 
    respTW = reshape(semrespV,9,8); semcomp2Act = respTW';
    
    mrespActL = NaN;  semrespActL = NaN; 
    xtl = {''}; trialsStr = cellfun(@num2str,trials','un',0);
    xaxL = NaN; xticks = []; xtickL =[];
    for ii = 1:8
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
        plot([ii ii],[4 35],'b-');
        if iii <= 8
            text(ii+2,35,sprintf('%s',xticklabels{iii}),'FontSize',6);
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
    xlim([0 length(mrespActL)+1]); ylim([3 45]);
    xlabel('Trial-Pairs');ylabel('Cells (%)');box off;
    set(gca,'xtick',xticks,'xticklabel',xtickL);
    legs = {'Conjunctive Cells      ','Complementary Cells 1','Complementary Cells 2',[9.5 0.1 43 0.2]}; 
    putLegendH(gca,legs,{'m','c','k'},'sigR',{[],'anova',[],6});
    format_axes(gca);
    changePosition(gca,[-0.08 0.1 0.17 -0.1]);
    save_pdf(hf,mData.pdf_folder,sprintf('tria_to_trial_unique_%c.pdf',G),600);
break;
end

%% 1 off diagnoal (uniqe between adjacent trials)
while 1
    %%
    respV = [];
    for an = 1:5
        respV(an,:) = diag(all_CI_mat(:,:,an));
    end

    
    mrespV = mean(respV); semrespV = std(respV)./sqrt(5);
    respTW = reshape(mrespV,10,8); mrespAct = respTW';
    respTW = reshape(semrespV,10,8); semrespAct = respTW';
    mrespActL = NaN;  semrespActL = NaN; 
    xtl = {''}; trialsStr = cellfun(@num2str,trials','un',0);
    xaxL = NaN; xticks = []; xtickL =[];
    for ii = 1:8
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
        plot([ii ii],[11 66],'b-');
        if iii <= 8
        text(ii+2,66,sprintf('%s',xticklabels{iii}),'FontSize',6);
        indsS = (theinds(iii)+1):(theinds(iii+1)-1);
        shadedErrorBar(indsS,mrespActL(indsS),semrespActL(indsS),{},0.5)
        iii=iii+1;
        end
        
    end
    xlim([0 length(mrespActL)+1]); ylim([10 70]);
    xlabel('Trials');ylabel('Cells (%)');box off;
    set(gca,'xtick',xticks,'xticklabel',xtickL);
%     legs = {'Responsive Cells',[9.5 0.1 34 0.2]}; 
%     putLegendH(gca,legs,{'k'},'sigR',{[],'anova',[],6});
    format_axes(gca);
    changePosition(gca,[-0.08 0.1 0.17 -0.1]);
    save_pdf(hf,mData.pdf_folder,sprintf('tria_to_trial_unique%c.pdf',G),600);
    %%
break;
end

%% anova running
while 1
    %%
    aVar = [];
    for an = 1:5
        tvar = conjV(an,:);
        tvar(10:10:79) = NaN; tvar(isnan(tvar)) = []; 
        aVar(an,:) = tvar;
    end
    [within,dvn,xlabels] = make_within_table({'cond','TrialsP'},[8,9]);
    dataT = make_between_table({aVar},dvn);
    rac = RMA(dataT,within);
    rac.ranova
    
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,rac,{'cond','hsd'},[1.5 1 1]);
%     h(h==1) = 0;
	xdata = make_xdata([8],[1 1.5]);
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
    save_pdf(hf,mData.pdf_folder,sprintf('conj_conds.pdf'),600);
    
    %%
    aVar = [];
    for an = 1:5
        tvar = comp1V(an,:);
        tvar(10:10:79) = NaN; tvar(isnan(tvar)) = []; 
        aVar(an,:) = tvar;
    end
    [within,dvn,xlabels] = make_within_table({'cond','TrialsP'},[8,9]);
    dataT = make_between_table({aVar},dvn);
    rac1 = RMA(dataT,within);
    rac1.ranova
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,rac1,{'cond','hsd'},[1.5 1 1]);
%     h(h==1) = 0;
	xdata = make_xdata([8],[1 1.5]);
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
    save_pdf(hf,mData.pdf_folder,sprintf('comp1_conds.pdf'),600);
    %%
    aVar = [];
    for an = 1:5
        tvar = comp2V(an,:);
        tvar(10:10:79) = NaN; tvar(isnan(tvar)) = []; 
        aVar(an,:) = tvar;
    end
    [within,dvn,xlabels] = make_within_table({'cond','TrialsP'},[8,9]);
    dataT = make_between_table({aVar},dvn);
    rac2 = RMA(dataT,within);
    rac2.ranova
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,rac2,{'cond','hsd'},[1.5 1 1]);
%     h(h==1) = 0;
	xdata = make_xdata([8],[1 1.5]);
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