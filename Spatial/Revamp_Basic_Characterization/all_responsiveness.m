function all_responsiveness
mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size;
ei = evalin('base','d15'); Rs = evalin('base','raster_data'); 
sC = evalin('base','selContexts'); rN = evalin('base','rasterNames');
an = 2;
cai_sampling_rate = ei{an}.thorExp.frameRate;
effective_sampling_rate = 1/0.15;
samplingRate = {'Ca','Ef','Ef','Ef','Ef','Ef','Ef','Ef','Ca','Ef'};
timeBefore = [2 5 7 NaN 7 NaN 7 NaN 2 5];

filename = fullfile(mData.pd_folder,'rasters.mat');
if 0
for an = 1:5
    anei = ei{an};
    anRs = Rs(an,:);
    for ii = 1:length(sC)
        sci = sC(ii);
        tRs = anRs{ii};
        thisContext = anei.plane{1}.contexts(sci);
        disp(thisContext.name)
        isCell = logical(tRs.iscell);
        if strcmp(samplingRate{ii},'Ca')
            tempR = tRs.fromFrames.sp_rasters;%(:,:,isCell);
            tempDur = tRs.fromFrames.duration;
            SR = cai_sampling_rate;
        else
            tempR = tRs.sp_rasters1;%(:,:,isCell);
            tempDur = tRs.duration1;
            SR = effective_sampling_rate;
        end
        rasters{an,ii} = find_resp_mdata(tempR,tempDur,timeBefore(ii),SR,thisContext.name,anei,tRs,isCell);
    end
end
save(filename,'rasters');
else
    load(filename);
end
n = 0;
%%
ssC = [1 2 3 5 7 9 10];
for an = 1:5
    resp = [];
    for ii = 1:length(ssC)
        tR = rasters{an,ssC(ii)};
        resp(:,ii) = tR.resp.p < 0.05;
    end
    all_resp{an} = resp;
end
%%
average_OI = nan(7,7);
all_OI_mat = [];
for an = 1:5
    pop_map = [];
    OI = [];
    for rr = 1:7
        for cc = 1:7
            rC1 = all_resp{an}(:,rr);
            rC2 = all_resp{an}(:,cc);
            pop_map(rr,cc) = sum(rC1&rC2)/length(rC1);
            OI(rr,cc) = sum(rC1&rC2)/(sum(rC1)+sum(rC2)-sum(rC1&rC2))
        end
    end
    all_pop_map{an} = pop_map;
    mask = triu(ones(size(OI)),1);mask(mask==0) = NaN
    all_OI{an} = OI.*mask;
    if an == 1
        average_OI = all_OI{an};
    else
        average_OI = average_OI + all_OI{an};
    end
    all_OI_mat(:,:,an) = all_OI{an};
end
average_OI = average_OI/5;

general_imagesc({1,[1 2 15 4]},all_OI')
figure(100);clf;imagesc(average_OI); colorbar;

%%
bins = [1:7];
ahi = nan(5,7);
for an = 1:5
    thismap = all_resp{an};
    hi = sum(thismap)/size(thismap,1);
    ahi(an,:) = hi;
end

dataT = array2table(ahi);dataT.Properties.VariableNames = {'C1','C2','C3','C4','C5','C6','C7',};
within = table([1 2 3 4 5 6 7]');within.Properties.VariableNames= {'Cs'};within.Cs = categorical(within.Cs);
ra = repeatedMeasuresAnova(dataT,within,0.05);
 mVar = ra.est_marginal_means.Mean;semVar = ra.est_marginal_means.Formula_StdErr;
combs = ra.mcs.combs; p = ra.mcs.p; h = ra.mcs.p < 0.05;


    xdata = [1:7]; 
    colors = mData.colors(1); colors = repmat(colors,1,8);
    hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.25 1],'color','w');
    hold on;
%     tcolors = {colors{1};colors{1};colors{2};colors{2};colors{3};colors{3};colors{4};colors{4}};
    tcolors = {colors{1};colors{2};colors{3};colors{4};colors{1};colors{2};colors{3};colors{4};};
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.1,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.001,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.7,'sigLinesStartYFactor',0.1);
    % plot([0.5 11],[-0.5 0.5],'linewidth',1.5)
%     maxY = maxY + 5 + 2;
    set(gca,'xlim',[0.25 xdata(end)+0.75],'ylim',[0 maxY],'FontSize',6,'FontWeight','Bold','TickDir','out');
    xticks = xdata(1:end); xticklabels = {'0','1','2','3','4','5','6','7'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
%     xtickangle(30)
    changePosition(gca,[0.1 0.03 -0.04 -0.05]);
    put_axes_labels(gca,{[],[0 0 0]},{'frac',[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('average fraction responses versus conditions '),600);

%%
bins = [0:1:7];
ahi = nan(5,8);
for an = 1:5
    thismap = all_resp{an};
    sumResp = sum(thismap,2);
    [hi,bi] = hist(sumResp,bins);
    hi = hi/sum(hi);
    ahi(an,:) = hi;
    all_sumResp{an} = sumResp;
end

[mVar,semVar] = findMeanAndStandardError(ahi);
    combs = []; p = []; h = [];
%     row = [1 2]; ii = ismember(combs,row,'rows'); p(ii) = mcGroup{1,6}; h(ii) = 1; 
    % row = [5 6]; ii = ismember(combs,row,'rows'); p(ii) = mcTI{2,6}; h(ii) = 1; 
    % row = [3 4]; ii = ismember(combs,row,'rows'); p(ii) = mcTI{3,6}; h(ii) = 1; 

    xdata = [1:8]; 
    colors = mData.colors(1); colors = repmat(colors,1,8);
    hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.25 1],'color','w');
    hold on;
%     tcolors = {colors{1};colors{1};colors{2};colors{2};colors{3};colors{3};colors{4};colors{4}};
    tcolors = {colors{1};colors{2};colors{3};colors{4};colors{1};colors{2};colors{3};colors{4};};
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',2,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.001,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.7,'sigLinesStartYFactor',0.1);
    % plot([0.5 11],[-0.5 0.5],'linewidth',1.5)
%     maxY = maxY + 5 + 2;
    set(gca,'xlim',[0.25 xdata(end)+0.75],'ylim',[0 maxY],'FontSize',6,'FontWeight','Bold','TickDir','out');
    xticks = xdata(1:end); xticklabels = {'0','1','2','3','4','5','6','7'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
%     xtickangle(30)
    changePosition(gca,[0.1 0.03 -0.04 -0.05]);
    put_axes_labels(gca,{[],[0 0 0]},{'frac',[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('average fraction responses versus conditions '),600);