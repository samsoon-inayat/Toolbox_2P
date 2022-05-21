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
        duration_onset_moveT = []; timeToCompleteT = [];
        for ii = 1:length(markers1)
            st = markers1(ii);
            se = markers2(ii);
            sp{ii} = speed(st:se);
            t{ii} = ts(st:se)-ts(st);
            ind(ii) = find((st:se)-markers2i(ii)>0,1,'first');
            t_on_move = find(sp{ii} > speed_threshold,1,'first');
            duration_onset_move(ii,cc,an) = t{ii}(t_on_move);
            duration_onset_moveT(ii) = t{ii}(t_on_move);
            
            timeToCompleteT(ii) = ts(markers2i(ii)) - ts(markers1i(ii));
            distT_cov(ii) = ds(markers2i(ii)) - ds(markers1i(ii));
            distI_cov(ii) = ds(markers2iI(ii)) - ds(markers1iI(ii));
        end
        duration_onset_moveC(an,cci) = mean(duration_onset_moveT);
        timeToCompleteC(an,cci) = mean(timeToCompleteT);
        timeToCompleteC_All(an,:) = (timeToCompleteT);
        distT(an,:) = (distT_cov);
        distI(an,:) = (distI_cov);
    end
end
n=0;
%%
if 0
moas = duration_onset_moveC;
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
% % xdata = [1 2 3]; maxY = 50;
% 
% xdata = [1:1.15:6]; xdata = xdata(1:3);
% maxY = 1;
colors = mData.colors(3:end);
hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.25 1],'color','w');
hold on;
tcolors = colors;
[hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
    'ySpacing',0.21,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.005,...
    'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.7,'sigLinesStartYFactor',0.051);
set(gca,'xlim',[0.25 max(xdata)+0.75],'ylim',[0 maxY],'FontSize',6,'FontWeight','Normal','TickDir','out');
xticks = xdata; xticklabels = {'Ar','ArL','Ar*'};
set(gca,'xtick',xticks,'xticklabels',xticklabels);
xtickangle(30);
changePosition(gca,[0.19 0.05 -0.4 -0.1])
put_axes_labels(gca,{'',[0 0 0]},{{'Movement','Latency (sec)'},[0 0 0]});
save_pdf(hf,mData.pdf_folder,'Duration to Movement Onset',600);
return
end
%%

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
xticks = xdata; xticklabels = {'Ar','ArL','Ar*'};
set(gca,'xtick',xticks,'xticklabels',xticklabels);
xtickangle(30);
changePosition(gca,[0.13 0.05 -0.4 -0.1])
put_axes_labels(gca,{'',[0 0 0]},{{'Trial Time (sec)'},[0 0 0]});
save_pdf(hf,mData.pdf_folder,'Time to Complete Trial',600);
