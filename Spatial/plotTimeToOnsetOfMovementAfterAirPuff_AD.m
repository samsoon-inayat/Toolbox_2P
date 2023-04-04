function plotTimeToOnsetOfMovementAfterAirPuff(b,markers1,markers2,fn)
%%
n = 0;
%%
ei_C = evalin('base','ei10_C');
ei_A = evalin('base','ei10_A');
mData = evalin('base','mData');

dt_C = get_duration_and_time(ei_C);
dt_A = get_duration_and_time(ei_A);
n = 0;
%%
if 1
moas = dt_C.duration;
moas_A = dt_A.duration;
for ii = 1:size(moas,2)
    varNames{ii} = sprintf('Trials_Cond%d',ii);
end
data = [moas;moas_A];
dataT = table([ones(size(moas,1),1);2*ones(size(moas_A,1),1)],data(:,1),data(:,2),data(:,3),data(:,4));
dataT.Properties.VariableNames = {'Group' varNames{1} varNames{2} varNames{3} varNames{4}};
dataT.Group = categorical(dataT.Group);
within = table([1 2 3 4]');
within.Properties.VariableNames = {'Condition'};
within.Condition = categorical(within.Condition);
ra = repeatedMeasuresAnova(dataT,within,0.05);

mVar = ra.est_marginal_means.Mean;
semVar = ra.est_marginal_means.Formula_StdErr;
combs = ra.mcs.combs; p = ra.mcs.p; h = ra.mcs.p < -0.05;

xdata = [1 2 3 4 6 7 8 9]; 
maxY = 3;
colors = mData.colors;
hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 2 1],'color','w');
hold on;
tcolors = colors(1:4);
tcolors = repmat(tcolors,2,1)';
hbs = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
    'maxY',maxY,'ySpacing',1,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
    'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.7,'sigLinesStartYFactor',-0.1);
set(gca,'xlim',[0.25 max(xdata)+0.75],'ylim',[0 maxY],'FontSize',6,'FontWeight','Bold','TickDir','out','xcolor','k','ycolor','k');
xticks = xdata; xticklabels = {'C1','C2','C3','C4'};
set(gca,'xtick',xticks,'xticklabels',xticklabels);
changePosition(gca,[0.17 0.02 -0.1 -0.011])
put_axes_labels(gca,{'Conditions',[0 0 0]},{{'Movement','Latency (s)'},[0 0 0]});
save_pdf(hf,mData.pdf_folder,'Duration to Movement Onset',600);
return
end
%%

moas = timeToCompleteC;
for ii = 1:size(moas,2)
    varNames{ii} = sprintf('Trials_Cond%d',ii);
end
data = [moas];
dataT = table(data(:,1),data(:,2),data(:,3),data(:,4));
dataT.Properties.VariableNames = {varNames{1} varNames{2} varNames{3} varNames{4}};
within = table([1 2 3 4]');
within.Properties.VariableNames = {'Condition'};
within.Condition = categorical(within.Condition);

% writetable(between,'Training_Data.xls');
rm = fitrm(dataT,'Trials_Cond1,Trials_Cond2,Trials_Cond3,Trials_Cond4~1','WithinDesign',within,'WithinModel','Condition');
rtable = ranova(rm,'WithinModel',rm.WithinModel);
mauchlytbl = mauchly(rm);
% multcompare(rm,'Day','ComparisonType','bonferroni')
mcTI = find_sig_mctbl(multcompare(rm,'Condition','ComparisonType','bonferroni'));

[mVar,semVar] = findMeanAndStandardError(moas);
combs = nchoosek(1:8,2); p = ones(size(combs,1),1); h = logical(zeros(size(combs,1),1));

xdata = [1:1.15:6]; xdata = xdata(1:4);
maxY = 25;
colors = mData.colors;
hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1 1],'color','w');
hold on;
tcolors = colors;
hbs = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
    'maxY',maxY,'ySpacing',1,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
    'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.7,'sigLinesStartYFactor',-0.1);
set(gca,'xlim',[0.25 max(xdata)+0.75],'ylim',[0 maxY],'FontSize',6,'FontWeight','Bold','TickDir','out');
xticks = xdata; xticklabels = {'C1','C2','C3','C4'};
set(gca,'xtick',xticks,'xticklabels',xticklabels);
changePosition(gca,[0.17 0.02 -0.1 -0.011])
put_axes_labels(gca,{'Conditions',[0 0 0]},{{'Trial Time (sec)'},[0 0 0]});
save_pdf(hf,mData.pdf_folder,'Time to Complete Trial',600);

function out = get_duration_and_time(ei)

timeBefore = 0;
timeAfter = 15;
speed_threshold = 1;
for an = 1:length(ei)
    for cc = 1:4
        b = ei{an}.b;
        markers1i = ei{an}.plane{1}.contexts(cc).markers.air_onsets;
        markers2i = ei{an}.plane{1}.contexts(cc).markers.air_offsets;
        fn = 101;
        speed = b.fSpeed;
        % figure(1000);clf;plot(b.speed);
        ts = b.ts;
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
        end
        duration_onset_moveC(an,cc) = mean(duration_onset_moveT);
        timeToCompleteC(an,cc) = mean(timeToCompleteT);
    end
end
out.duration = duration_onset_moveC;
out.timeToComplete = timeToCompleteC;