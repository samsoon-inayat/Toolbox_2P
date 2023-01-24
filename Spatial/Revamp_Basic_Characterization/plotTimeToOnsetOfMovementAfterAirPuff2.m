function plotTimeToOnsetOfMovementAfterAirPuff(b,markers1,markers2,fn)
%%
n = 0;
%%
ei1 = evalin('base','ei');
ei = [ei1(1:5)];
speed_threshold = 0;

mData = evalin('base','mData');

timeBefore = 0;
timeAfter = 15;
cs = [3 4 5];
for an = 1:length(ei)
    for cci = 1:3
        cc = cs(cci);
        b = ei{an}.b;
        if cci == 2
            markers1i_L = ei{an}.plane{1}.contexts(4).markers.light22_onsets;
        end
        markers1i = ei{an}.plane{1}.contexts(cc).markers.air_onsets;
        markers2i = ei{an}.plane{1}.contexts(cc).markers.air_offsets;
        markers1iI = markers2i(1:end);
        markers2iI = markers1i(2:end);
        ts = b.ts;
        meandur = mean(ts(markers2iI) - ts(markers1iI(1:9)));
        markers2iI(10) = markers1iI(end) + round(1e6 * meandur/b.si);
        fn = 101;
        speed = b.fSpeed;
        % figure(1000);clf;plot(b.speed);
        
        ds = b.dist;
        markers1 = markers1i - round(1e6 * timeBefore/b.si);
        markers2 = markers2i + round(1e6 * timeAfter/b.si);
        markers3 = markers2i + round(1e6 * 7.5/b.si);
        markers1m1 = markers1i - round(1e6 * 0/b.si);
        markers2m1 = markers2i + round(1e6 * 0/b.si);
        duration_onset_moveT = []; timeToCompleteT = [];
        if cci == 2
            markers1m1_L = markers1i_L - round(1e6 * 0/b.si);
        end
        
        for ii = 1:length(markers1)
            st = markers1(ii);
            se = markers2(ii);
            sp{ii} = speed(st:se);
            t{ii} = ts(st:se)-ts(st);
            ind(ii) = find((st:se)-markers2i(ii)>0,1,'first');
            t_on_move = find(sp{ii} > speed_threshold,1,'first');
            duration_onset_move(ii,cci,an) = t{ii}(t_on_move);
            duration_onset_moveT(ii) = t{ii}(t_on_move);
            speed_at_onset(ii,cci,an) = speed(markers1m1(ii));
            speed_at_offset(ii,cci,an) = speed(markers2m1(ii));
            speed_at_offset_c(ii,cci,an) = speed(markers3(ii));
            if cci == 2
                speed_at_onset_L(ii,an) = speed(markers1m1_L(ii));
            end
            
            timeToCompleteT(ii) = ts(markers2i(ii)) - ts(markers1i(ii));
            distT_cov(ii) = ds(markers2i(ii)) - ds(markers1i(ii));
            distI_cov(ii) = ds(markers2iI(ii)) - ds(markers1iI(ii));
            distI_cov_all(ii,cci,an) = ds(markers2iI(ii)) - ds(markers1iI(ii));
        end
        duration_onset_moveC(an,cci) = mean(duration_onset_moveT);
        timeToCompleteC(an,cci) = mean(timeToCompleteT);
        timeToCompleteC_All(an,:) = (timeToCompleteT);
        distT(an,:) = (distT_cov);
        distI(an,:) = (distI_cov);
    end
end
dom = reshape(duration_onset_move,30,5); dom1 = dom; dom1(dom==0) = nan; 
out = descriptiveStatistics(dom1(:),'decimal_places',2); 
outpim = descriptiveStatistics(100*(sum(isnan(dom1))'/30),'decimal_places',2);
sao = reshape(speed_at_onset,30,5);
saof = reshape(speed_at_offset,30,5);
saof_c = reshape(speed_at_offset_c,30,5);
dI = reshape(distI_cov_all,30,5);

sao_L = speed_at_onset_L;

sao_L1 = mean(speed_at_onset_L,1);

sao1 = squeeze(mean(speed_at_onset,1));
saof1 = squeeze(mean(speed_at_offset,1));
saof1_c = squeeze(mean(speed_at_offset_c,1));

sase = [sao1' saof1' saof1_c'];
[within,dvn,xlabels] = make_within_table({'SE','Cond'},[3,3]);
dataT = make_between_table({sase},dvn);
ra = RMA(dataT,within,{'tukey-kramer','hsd'});
ra.ranova
n=0;

%%

out = descriptiveStatistics(speed_at_onset_L(:),'decimal_places',2);
out = descriptiveStatistics(speed_at_onset(:),'decimal_places',2);
out = descriptiveStatistics(speed_at_offset(:),'decimal_places',2);
out = descriptiveStatistics(speed_at_offset_c(:),'decimal_places',2);
%% movemenet latency distrubtion
while 1
    clear distD
    for an = 1:5
        distD{an,1} = dom(:,an);
    end
    
    [distDo,allVals] = getAveragesAndAllValues(distD);
   minBin = min(allVals);
   maxBin = max(allVals);
   incr = 0.25;
   hf = figure(8);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.5 1],'color','w');
   hold on;
   [ha,hb,hca] = plotAverageDistributions(distD,'colors',{'k'},'maxY',100,'min',minBin,'incr',incr,'max',maxBin,'pdf_or_cdf','pdf');
   set(gca,'FontSize',6,'FontWeight','Normal','TickDir','out','xcolor','k','ycolor','k');
   xlim([minBin maxBin]);
   changePosition(gca,[0.12 0.13 -0.2 -0.13]);
    put_axes_labels(gca,{'Latency (s)',[0 0 0]},{{'Trials (%)'},[0 0 0]});
    format_axes(gca);
    save_pdf(hf,mData.pdf_folder,sprintf('Distribution_movement_latency'),600);
    %%
    break;
end

%% speed at air onset or offset or light onset (change variable to do so)
while 1
    ff = makeFigureRowsCols(107,[1 0.5 4 0.5],'RowsCols',[1 4],'spaceRowsCols',[0.1 0.07],'rightUpShifts',[0.08 0.3],'widthHeightAdjustment',[-80 -400]);
    set(gcf,'color','w'); set(gcf,'Position',[5 5 4 1]);
    vars = {'sao','saof','sao_L','saof_c'};
    for ii = 1
    clear distD
    for an = 1:5
        cmdTxt = sprintf('distD{an,1} = %s(:,an);',vars{ii}); eval(cmdTxt);
        avgDist(an) = 100*(sum(distD{an,1}>0)/30);
    end
    
    [distDo,allVals] = getAveragesAndAllValues(distD);
   minBin = min(allVals);
   maxBin = max(allVals);
   incr = 3;
%    hf = figure(8);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.1 1],'color','w');
    axes(ff.h_axes(1,ii));
   hold on;
   [ha,hb,hca] = plotAverageDistributions(distD,'colors',{'k'},'maxY',100,'min',minBin,'incr',incr,'max',maxBin,'pdf_or_cdf','cdf');
   set(gca,'FontSize',6,'FontWeight','Normal','TickDir','out','xcolor','k','ycolor','k');
   xlim([0 maxBin]); ylim([0 100]);
%    changePosition(gca,[0.12 0.13 -0.0 -0.13]);
   set(gca,'xtick',[0 10 20 30])
   if ii == 1
        put_axes_labels(gca,{'Speed (cm/s)',[0 0 0]},{{'Trials (%)'},[0 0 0]});
   else
       put_axes_labels(gca,{'Speed (cm/s)',[0 0 0]},{{''},[0 0 0]});
   end
    format_axes(gca);
    end
    save_pdf(gcf,mData.pdf_folder,sprintf('Distribution'),600);
    %%
    break;
end

%% distribution distance covered
while 1
    clear distD
    for an = 1:5
        distD{an,1} = dI(:,an);
    end
    
    [distDo,allVals] = getAveragesAndAllValues(distD);
   minBin = min(allVals);
   maxBin = max(allVals);
   incr = 10;
   hf = figure(8);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.5 1],'color','w');
   hold on;
   [ha,hb,hca] = plotAverageDistributions(distD,'colors',{'k'},'maxY',100,'min',minBin,'incr',incr,'max',maxBin,'pdf_or_cdf','cdf');
   set(gca,'FontSize',6,'FontWeight','Normal','TickDir','out','xcolor','k','ycolor','k');
   xlim([minBin maxBin]);
   changePosition(gca,[0.12 0.13 -0.2 -0.13]);
    put_axes_labels(gca,{'Distance (cm)',[0 0 0]},{{'Trials (%)'},[0 0 0]});
    format_axes(gca);
    save_pdf(hf,mData.pdf_folder,sprintf('Distribution'),600);
    %%
    break;
end


%% average speed at onset versus offset RMANOVA
while 1
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'SE','hsd'},[1.5 1 1]);
    xdata = make_xdata([3 3 3],[1 1.5]);
    hf = get_figure(5,[8 7 1.5 1]);
    tcolors = repmat(mData.colors(1:3),1,3);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',5,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.15);
    maxY = maxY;
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1)-1 maxY+1]); format_axes(gca);
    xticks = xdata; xticklabels = {'C3','C4','C5'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0.13 0.01 -0.4 0]); put_axes_labels(gca,{[],[0 0 0]},{{'Average Speed','(cm/s)'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('speed.pdf'),600);
break;
end
%%
while 1
dataD = {dom(:,1);dom(:,2);dom(:,3);dom(:,4);dom(:,5)};
[distDo,allVals] = getAveragesAndAllValues(dataD);
   minBin = min(allVals);
   maxBin = max(allVals);
   incr = 0.25;
   hf = figure(8);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.25 1],'color','w');
   hold on;
   [ha,hb,hca] = plotAverageDistributions(dataD,'colors',{'k'},'maxY',100,'min',minBin,'incr',incr,'max',maxBin,'pdf_or_cdf','pdf');
   set(gca,'FontSize',6,'FontWeight','Normal','TickDir','out','xcolor','k','ycolor','k');
    xlabel('Movement Latency (s)');
    ylabel('%');
    changePosition(gca,[0.01 0.01 0 0]);
    save_pdf(hf,mData.pdf_folder,'movement latency distribution',600);
break;
end
%% anova movement latency across conditions
while 1
moas = duration_onset_moveC;
moas = squeeze(mean(duration_onset_move,1));
moas11 = reshape(duration_onset_move,30,5); moas11(moas11 == 0) = NaN;
moas = moas';
for ii = 1:size(moas,2)
    varNames{ii} = sprintf('Trials_Cond%d',ii);
end
data = [moas];
[within,dvn,xlabels] = make_within_table({'Cond'},[3]);
dataT = make_between_table({data},dvn);
ra = RMA(dataT,within,{'tukey-kramer','hsd'});
ra.ranova
[xdata,mVar,semVar,combs,p,h,colorsi,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Cond','hsd'},[1 1 1]);

colors = mData.colors(3:end);
hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.25 1],'color','w');
hold on;
tcolors = colors;
[hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
    'ySpacing',0.21,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.005,...
    'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.7,'sigLinesStartYFactor',0.051);
set(gca,'xlim',[0.25 max(xdata)+0.75],'ylim',[0 maxY],'FontSize',6,'FontWeight','Normal','TickDir','out');
xticks = xdata; xticklabels = {'C3','C4','C5'};
set(gca,'xtick',xticks,'xticklabels',xticklabels);
xtickangle(30);
changePosition(gca,[0.19 0.05 -0.4 -0.1])
put_axes_labels(gca,{'',[0 0 0]},{{'Movement','Latency (sec)'},[0 0 0]});
save_pdf(hf,mData.pdf_folder,'Duration to Movement Onset',600);
break;
end

%% anova distance covered during intertrials across conditions
while 1
moas = squeeze(mean(distI_cov_all,1));
moas = moas';
for ii = 1:size(moas,2)
    varNames{ii} = sprintf('Trials_Cond%d',ii);
end
data = [moas];
[within,dvn,xlabels] = make_within_table({'Cond'},[3]);
dataT = make_between_table({data},dvn);
ra = RMA(dataT,within,{'tukey-kramer','hsd'});
ra.ranova
[xdata,mVar,semVar,combs,p,h,colorsi,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Cond','hsd'},[1 1 1]);

colors = mData.colors(3:end);
hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.25 1],'color','w');
hold on;
tcolors = colors;
[hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
    'ySpacing',0.21,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.005,...
    'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.7,'sigLinesStartYFactor',0.051);
set(gca,'xlim',[0.25 max(xdata)+0.75],'ylim',[0 maxY],'FontSize',6,'FontWeight','Normal','TickDir','out');
xticks = xdata; xticklabels = {'C3','C4','C5'};
set(gca,'xtick',xticks,'xticklabels',xticklabels);
xtickangle(30);
changePosition(gca,[0.19 0.05 -0.4 -0.1])
put_axes_labels(gca,{'',[0 0 0]},{{'Distance','Traveled (cm)'},[0 0 0]});
save_pdf(hf,mData.pdf_folder,'Duration to Movement Onset',600);
break;
end
%% Trial Time
while 1
    moas = timeToCompleteC;
    for ii = 1:size(moas,2)
        varNames{ii} = sprintf('Trials_Cond%d',ii);
    end
    data = [moas];
    [within,dvn,xlabels] = make_within_table({'Cond'},[3]);
    dataT = make_between_table({data},dvn);
    ra = RMA(dataT,within,{'tukey-kramer','hsd'});
    ra.ranova
    [xdata,mVar,semVar,combs,p,h,colorsi,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Cond','hsd'},[1 1 1]);
    % dataT = array2table(data);
    % dataT.Properties.VariableNames = {varNames{1} varNames{2} varNames{3}};
    % within = table([1 2 3]');
    % within.Properties.VariableNames = {'Condition'};
    % within.Condition = categorical(within.Condition);
    % 
    % ra = repeatedMeasuresAnova(dataT,within);
    % 
    % mVar = ra.est_marginal_means.Mean;
    % semVar = ra.est_marginal_means.Formula_StdErr;
    % combs = ra.mcs.combs; p = ra.mcs.p; h = ra.mcs.p < 0.05;
    % 
    % xdata = [1:1.15:6]; xdata = xdata(1:3);
    % maxY = 15;
    colors = mData.colors(3:end);
    hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.25 1],'color','w');
    hold on;
    tcolors = colors;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',1,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.7,'sigLinesStartYFactor',-0.1);
    set(gca,'xlim',[0.25 max(xdata)+0.75],'ylim',[0 maxY],'FontSize',6,'FontWeight','Normal','TickDir','out');
    xticks = xdata; xticklabels = {'C3','C4','C5'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
    xtickangle(30);
    changePosition(gca,[0.13 0.05 -0.4 -0.1])
    put_axes_labels(gca,{'',[0 0 0]},{{'Trial Time (sec)'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,'Time to Complete Trial',600);
%%
break;
end
