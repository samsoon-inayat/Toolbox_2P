function trial_to_trial_Analysis_brake_vs_nobrake
%%
si = [Lb Ab_On Ab_Off Ar_On Ar_Off ArL_On ArL_Off Ars_On Ars_Off Lbs Abs_On Abs_Off ArL_L];
all_resp_A = find_all_trials_resp(oA,si);
all_resp_C = find_all_trials_resp(oC,si);
trials = mat2cell([1:10]',ones(size([1:10]')));
%%
[respRVA,conjVA,comp1VA,comp2VA] = find_resp_conj_comp(all_resp_A);
[respRVC,conjVC,comp1VC,comp2VC] = find_resp_conj_comp(all_resp_C);

%%
    aVar = []; aVarA = [];
    for an = 1:5
        tvar = conjVC(an,:);
        tvar(10:10:129) = NaN; tvar(isnan(tvar)) = []; 
        aVar(an,:) = tvar;
        tvar = conjVA(an,:);
        tvar(10:10:129) = NaN; tvar(isnan(tvar)) = []; 
        aVarA(an,:) = tvar;
    end
    tvar = conjVC(6,:);
    tvar(10:10:129) = NaN; tvar(isnan(tvar)) = []; 
    aVar(6,:) = tvar;
    [within,dvn,xlabels] = make_within_table({'cond','TrialsP'},[13,9]);
    dataT = make_between_table({aVar;aVarA},dvn);
    rac = RMA(dataT,within);
    rac.ranova
  
 %%
    aVar = []; aVarA = [];
    for an = 1:5
        tvar = comp1VC(an,:);
        tvar(10:10:129) = NaN; tvar(isnan(tvar)) = []; 
        aVar(an,:) = tvar;
        tvar = comp1VA(an,:);
        tvar(10:10:129) = NaN; tvar(isnan(tvar)) = []; 
        aVarA(an,:) = tvar;
    end
    tvar = comp1VC(6,:);
    tvar(10:10:129) = NaN; tvar(isnan(tvar)) = []; 
    aVar(6,:) = tvar;
    [within,dvn,xlabels] = make_within_table({'cond','TrialsP'},[13,9]);
    dataT = make_between_table({aVar;aVarA},dvn);
    rac1 = RMA(dataT,within);
    rac1.ranova
    
 %%
    aVar = []; aVarA = [];
    for an = 1:5
        tvar = comp2VC(an,:);
        tvar(10:10:129) = NaN; tvar(isnan(tvar)) = []; 
        aVar(an,:) = tvar;
        tvar = comp2VA(an,:);
        tvar(10:10:129) = NaN; tvar(isnan(tvar)) = []; 
        aVarA(an,:) = tvar;
    end
    tvar = comp2VC(6,:);
    tvar(10:10:129) = NaN; tvar(isnan(tvar)) = []; 
    aVar(6,:) = tvar;
    [within,dvn,xlabels] = make_within_table({'cond','TrialsP'},[13,9]);
    dataT = make_between_table({aVar;aVarA},dvn);
    rac2 = RMA(dataT,within);
    rac2.ranova

%% 1 off diagnoal (uniqe between adjacent trials)
while 1
    xticklabels = rasterNamesTxt(si);
    
    respRV = respRVA; conjV = conjVA; comp1V = comp1VA; comp2V = comp2VA;
%     respRV = respRVC; conjV = conjVC; comp1V = comp1VC; comp2V = comp2VC;
    
    respV = respRV;
    mrespV = mean(respV); semrespV = std(respV)./sqrt(5);
    respTW = reshape(mrespV,10,13); mrespAct = respTW'; 
    respTW = reshape(semrespV,10,13); semrespAct = respTW';
    
    respV = conjV;
    mrespV = mean(respV); semrespV = std(respV)./sqrt(5);
    mrespV(10:10:129) = NaN; mrespV(130) = NaN; mrespV(isnan(mrespV)) = []; 
    semrespV(10:10:129) = NaN; semrespV(130) = NaN; semrespV(isnan(semrespV)) = []; 
    respTW = reshape(mrespV,9,13); mconjAct = respTW'; 
    respTW = reshape(semrespV,9,13); semconjAct = respTW';
    
    respV = comp1V;
    mrespV = mean(respV); semrespV = std(respV)./sqrt(5);
    mrespV(10:10:129) = NaN; mrespV(130) = NaN; mrespV(isnan(mrespV)) = []; 
    semrespV(10:10:129) = NaN; semrespV(130) = NaN; semrespV(isnan(semrespV)) = []; 
    respTW = reshape(mrespV,9,13); mcomp1Act = respTW'; 
    respTW = reshape(semrespV,9,13); semcomp1Act = respTW';
    
    respV = comp2V;
    mrespV = mean(respV); semrespV = std(respV)./sqrt(5);
    mrespV(10:10:129) = NaN; mrespV(130) = NaN; mrespV(isnan(mrespV)) = []; 
    semrespV(10:10:129) = NaN; semrespV(130) = NaN; semrespV(isnan(semrespV)) = []; 
    respTW = reshape(mrespV,9,13); mcomp2Act = respTW'; 
    respTW = reshape(semrespV,9,13); semcomp2Act = respTW';
    
    mrespActL = NaN;  semrespActL = NaN; 
    xtl = {''}; trialsStr = cellfun(@num2str,trials','un',0);
    xaxL = NaN; xticks = []; xtickL =[];
    for ii = 1:13
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
        if iii <= 13
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

%% 1 off diagnoal (uniqe between adjacent trials)
while 1
    %%
    respV = [];
    for an = 1:5
        respV(an,:) = diag(all_CI_mat(:,:,an));
    end
    
    [within,dvn,xlabels] = make_within_table({'cond','Trials'},[13,10]);
    dataT = make_between_table({respV},dvn);
    ra = RMA(dataT,within);
    ra.ranova
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'cond','hsd'},[1.5 1 1]);
%     h(h==1) = 0;
	xdata = make_xdata([13],[1 1.5]);
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
    respTW = reshape(mrespV,10,13); mrespAct = respTW';
    respTW = reshape(semrespV,10,13); semrespAct = respTW';
    mrespActL = NaN;  semrespActL = NaN; 
    xtl = {''}; trialsStr = cellfun(@num2str,trials','un',0);
    xaxL = NaN; xticks = []; xtickL =[];
    for ii = 1:13
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
        if iii <= 13
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
        tvar(10:10:129) = NaN; tvar(isnan(tvar)) = []; 
        aVar(an,:) = tvar;
    end
    [within,dvn,xlabels] = make_within_table({'cond','TrialsP'},[13,9]);
    dataT = make_between_table({aVar},dvn);
    rac = RMA(dataT,within);
    rac.ranova
    
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,rac,{'cond','hsd'},[1.5 1 1]);
%     h(h==1) = 0;
	xdata = make_xdata([13],[1 1.5]);
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
        tvar(10:10:129) = NaN; tvar(isnan(tvar)) = []; 
        aVar(an,:) = tvar;
    end
    [within,dvn,xlabels] = make_within_table({'cond','TrialsP'},[13,9]);
    dataT = make_between_table({aVar},dvn);
    rac1 = RMA(dataT,within);
    rac1.ranova
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,rac1,{'cond','hsd'},[1.5 1 1]);
%     h(h==1) = 0;
	xdata = make_xdata([13],[1 1.5]);
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
        tvar(10:10:129) = NaN; tvar(isnan(tvar)) = []; 
        aVar(an,:) = tvar;
    end
    [within,dvn,xlabels] = make_within_table({'cond','TrialsP'},[13,9]);
    dataT = make_between_table({aVar},dvn);
    rac2 = RMA(dataT,within);
    rac2.ranova
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,rac2,{'cond','hsd'},[1.5 1 1]);
%     h(h==1) = 0;
	xdata = make_xdata([13],[1 1.5]);
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