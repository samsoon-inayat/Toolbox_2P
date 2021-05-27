function figure1_Light_Responsive

mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size;
ei = evalin('base','d15'); 
% selContexts = [1 2 3 4 4 5 6];
% rasterNames = {'light22T','airD','light22T','airD','airD','light22T'};
% Rs = get_rasters_data(ei,selContexts,rasterNames);

selContexts = [1 4 6 2 7 3 4 5];
rasterNames = {'light22T','light22T','light22T','air55T','air55T','airD','airD','airD'};
Rs = get_rasters_data(ei,selContexts,rasterNames);

an = 1;
cai_sampling_rate = ei{an}.thorExp.frameRate;
effective_sampling_rate = 1/0.15;
samplingRate = {'Ca','Ca','Ef','Ca'};
timeBefore = [2 2 NaN 2];

mRs = calc_mean_rasters(Rs,1:10);
Rs = find_responsive_rasters(Rs,1:10);
view_population_vector(Rs,mRs,100);


for rr = 1:size(Rs,1)
    ccs = [];
    for cc = 1:size(Rs,2)
        R = Rs{rr,cc};
        ccs(:,cc) = R.resp.vals';
        resp_fraction(rr,cc) = R.resp.fraction;
    end
    resp_vals{rr} = ccs;
end

for ii = 1:length(resp_vals)
    ccs = resp_vals{ii};
    OI = NaN(size(ccs,2),size(ccs,2));
    mask = triu(ones(size(ccs,2)),1);
    for rr = 1:size(ccs,2)
        ccs1 = ccs(:,rr);
        for cc = 1:size(ccs,2)
            if mask(rr,cc)
                ccs2 = ccs(:,cc);
                shared = ccs1 & ccs2;
                OI(rr,cc) = sum(shared)/(sum(ccs1)+sum(ccs2)-sum(shared));
            end
        end
    end
    all_OI{ii} = OI;
end

figure(200);clf
for ii = 1:length(resp_vals)
    OI = all_OI{ii};
    subplot(5,1,ii);
    imagesc(OI);
    colorbar;
end


dataT = array2table(resp_fraction(:,1:3));
dataT.Properties.VariableNames = {'L1','L2','L3'};
within = array2table([1 2 3]');
within.Properties.VariableNames = {'Cond'};
within.Cond = categorical(within.Cond);
ra = repeatedMeasuresAnova(dataT,within,0.05);

%%
mVar = ra.est_marginal_means.Mean;semVar = ra.est_marginal_means.Formula_StdErr;
    combs = ra.mcs.combs; p = ra.mcs.p; h = ra.mcs.p < 0.05;
%     row = [1 2]; ii = ismember(combs,row,'rows'); p(ii) = mcGroup{1,6}; h(ii) = 1; 
    % row = [5 6]; ii = ismember(combs,row,'rows'); p(ii) = mcTI{2,6}; h(ii) = 1; 
    % row = [3 4]; ii = ismember(combs,row,'rows'); p(ii) = mcTI{3,6}; h(ii) = 1; 

    xdata = [1 2 3]; 
    colors = mData.colors;
    hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.25 1],'color','w');
    hold on;
%     tcolors = {colors{1};colors{1};colors{2};colors{2};colors{3};colors{3};colors{4};colors{4}};
    tcolors = {colors{1};colors{2};colors{3};colors{4};colors{1};colors{2};colors{3};colors{4};};
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.1,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.7,'sigLinesStartYFactor',0.1);
    for ii = 5:length(hbs)
    set(hbs(ii),'facecolor','none','edgecolor',tcolors{ii});
    end
    % plot([0.5 11],[-0.5 0.5],'linewidth',1.5)
%     maxY = maxY + 5 + 2;
    set(gca,'xlim',[0.25 xdata(end)+0.75],'ylim',[0 maxY],'FontSize',6,'FontWeight','Bold','TickDir','out');
    xticks = xdata(1:end); xticklabels = {'C1','C4','C1'''};
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
%     xtickangle(30)
    changePosition(gca,[0.08 0.03 -0.04 -0.05]);
    put_axes_labels(gca,{[],[0 0 0]},{'Fraction',[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('fraction_light_responsive'),600);
%%

return;
sci = 1;
for anii = 1:5
    for sci = 1:2
        sc = selContexts(sci);
        anei = ei{anii};
        tRs = Rs{anii,sci};
        tempR = tRs.fromFrames.sp_rasters;%(:,:,isCell);
        tempDur = tRs.fromFrames.duration;
        isCell = logical(tRs.iscell);
        thisContext = anei.plane{1}.contexts(selContexts(sci));
        SR = cai_sampling_rate;
        rasters{anii,sci} = find_resp_mdata(tempR,tempDur,timeBefore(sci),SR,thisContext.name,anei,tRs,isCell);
    end
end

for anii = 1:5
    for sci = 3
        sc = selContexts(sci);
        anei = ei{anii};
        tRs = Rs{anii,sci};
        tempR = tRs.sp_rasters1;%(:,:,isCell);
        tempDur = tRs.duration1;
        SR = effective_sampling_rate;
        isCell = logical(tRs.iscell);
        thisContext = anei.plane{1}.contexts(selContexts(sci));
        rastersD{anii,1} = find_resp_mdata(tempR,tempDur,NaN,SR,thisContext.name,anei,tRs,isCell);
    end
end
meanRs = calc_mean_rasters(rasters,1:10);
meanRs = correct_for_size(meanRs);
[all_ccs,fR] = get_responsive_cells(rasters);
[all_ccs,fR] = get_responsive_cells_zMI(rastersD);
meanRsD = calc_mean_rasters(rastersD,1:10);
%%
cell_list = all_ccs{an,2};
for an = 1:5
%     [popVs{an},cpc{an},ccc{an}] = calc_pop_vector_corr(meanRs(an,:),all_ccs{an,1} & all_ccs{an,2},[1 1]);
%     [popVsD{an},~,~,cell_nums] = calc_pop_vector_corr(meanRsD(an,:),cell_list,[]);
    [popVs{an},cpc{an},ccc{an},~] = calc_pop_vector_corr(meanRs(an,:),cell_list,[1 2]);
    
end
graphs = [];
for an = 1:5
    for cc = 1:2
        graphs{an,cc} = popVs{an}{cc}.popV;
    end
%     graphs{an,3} = popVsD{an}{1}.popV;
end

general_imagesc({1,[1 2 15,8]},graphs)
    
n = 0;
%% Average correlation
figure(1);clf;set(gcf,'Units','Inches');set(gcf,'Position',[1 2 7 3])
figure(2);clf;set(gcf,'Units','Inches');set(gcf,'Position',[1 2 7 3])
for anii = 1:5
%     figure(1);subplot(1,5,anii)
%     thisCorrV = popVs{anii}{1}.corrV; imagesc(thisCorrV);
%     figure(2);subplot(1,5,anii)
%     thisCorrV = popVs{anii}{1}.popV; imagesc(thisCorrV);
    
%     figure(1);subplot(1,5,anii)
%     thisCorrV = popVs{anii}{1}.popV; imagesc(thisCorrV);
%     figure(2);subplot(1,5,anii)
%     thisCorrV = popVs{anii}{2}.popV; imagesc(thisCorrV);
    figure(1);subplot(1,5,anii)
    thisCorrV = popVs{anii}{1}.popV; imagesc(thisCorrV);
    figure(2);subplot(1,5,anii)
    thisCorrV = popVs{anii}{2}.popV; imagesc(thisCorrV);

end

%%
runthis = 1;
if runthis
    dataT = array2table(fR);
    dataT.Properties.VariableNames = {'C1','C4'};
    within = table([1 2]');
    within.Properties.VariableNames = {'Condition'};
    within.Condition = categorical(within.Condition);
    ra = repeatedMeasuresAnova(dataT,within);
%%
    mVar = ra.est_marginal_means.Mean;semVar = ra.est_marginal_means.Formula_StdErr;
    combs = ra.mcs.combs; p = ra.mcs.p; h = ra.mcs.p < 0.05;
%     row = [1 2]; ii = ismember(combs,row,'rows'); p(ii) = mcGroup{1,6}; h(ii) = 1; 
    % row = [5 6]; ii = ismember(combs,row,'rows'); p(ii) = mcTI{2,6}; h(ii) = 1; 
    % row = [3 4]; ii = ismember(combs,row,'rows'); p(ii) = mcTI{3,6}; h(ii) = 1; 

    xdata = [1 2]; 
    colors = mData.colors;
    hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.25 1],'color','w');
    hold on;
%     tcolors = {colors{1};colors{1};colors{2};colors{2};colors{3};colors{3};colors{4};colors{4}};
    tcolors = {colors{1};colors{2};colors{3};colors{4};colors{1};colors{2};colors{3};colors{4};};
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.1,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.7,'sigLinesStartYFactor',0.1);
    for ii = 5:length(hbs)
    set(hbs(ii),'facecolor','none','edgecolor',tcolors{ii});
    end
    % plot([0.5 11],[-0.5 0.5],'linewidth',1.5)
%     maxY = maxY + 5 + 2;
    set(gca,'xlim',[0.25 xdata(end)+0.75],'ylim',[0 maxY],'FontSize',6,'FontWeight','Bold','TickDir','out');
    xticks = xdata(1:end); xticklabels = {'C1','C2'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
%     xtickangle(30)
    changePosition(gca,[0.02 0.03 -0.04 -0.05]);
    put_axes_labels(gca,{[],[0 0 0]},{'zMI',[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('light_responsive'),600);
return;
end

