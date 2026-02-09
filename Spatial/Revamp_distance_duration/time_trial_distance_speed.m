function time_trial_distance_speed

udata = evalin('base','udata1');
ei = evalin('base','ei');
configurations = {'C3','C4','C5'};
air_phases = {'ON','OFF'};

bin_types = {'time_bin','dist_bin'};
variable_combs = {'time_dist','time_speed','dist_speed','FR_time','FR_dist','FR_speed'};
variable_combs = {'time_dist','time_speed','dist_speed'};
ow = [2 2] ;

avar = [];
for an = 1:5
    data_an = udata{an};
    field_names = fieldnames(data_an);
    for ii = 1:length(field_names)
        varname = field_names{ii};
        cmdTxt = sprintf('%s = data_an.%s;',varname,varname);eval(cmdTxt);
    end
    anvar = [];
    for cn = 1:length(configurations)
        for ap = 1:2
            for tn = 1:10
                % get time to complete trial
                if ap == 1
                    cmdTxt = sprintf('air_trials = air_trials_on .* %s;',configurations{cn}); eval(cmdTxt);
                else
                    cmdTxt = sprintf('air_trials = air_trials_off .* %s;',configurations{cn}); eval(cmdTxt);
                end
                tts = ts(air_trials == tn);  tds = ds(air_trials == tn); tsp = speed(air_trials == tn);
                % select which variable we want to accumulate and run
                % RM-ANOVA - it could be time, distance, or speed. Also
                % choose ap 1 or 2 above in the for loop for choosing
                % either air-on or off phases or both
                
                % thisvar = tts(end) - tts(1);
                thisvar = tds(end) - tds(1);
                % thisvar = skewness(speed(air_trials == tn));
                % vn = 2; bn = 1;
                % params = {'no_of_bins_for_MI',10,'no_of_shuffles_for_norm',0,'animal_info',ei{an},'overwrite_processing',ow,'air_phase',air_phases{ap},...
                %         'configuration',configurations{cn},'variables',variable_combs{vn},'bin_type',bin_types{bn},'trial_type','individual'};
                %     met = myMetrics(tts',tsp',params);
                % 
                % temp0 = corr(tts',tds'); temp1 = corr(tts',tsp'); temp2 = corr(tds',tsp');
                % thisvar = met.PC(1);
                % anvar = [anvar thisvar];% outD.trial_metrics(:,idx)'];
                % 
                %  vn = 3; bn = 1;
                % params = {'no_of_bins_for_MI',10,'no_of_shuffles_for_norm',0,'animal_info',ei{an},'overwrite_processing',ow,'air_phase',air_phases{ap},...
                %         'configuration',configurations{cn},'variables',variable_combs{vn},'bin_type',bin_types{bn},'trial_type','individual'};
                %     met = myMetrics(tds',tsp',params);
                % 
                % % temp0 = corr(tts',tds'); temp1 = corr(tts',tsp'); temp2 = corr(tds',tsp');
                % thisvar = met.PC(1);
                anvar = [anvar thisvar];% outD.trial_metrics(:,idx)'];
                n = 0;
            end
        end
    end
    avar(an,:) = anvar;
end
n = 0;
%%
n = 0;
[within,dvn,xlabels,awithinD] = make_within_table({'CN','AP','TN'},[3,2,10]);
dataT = make_between_table({avar},dvn);
ra = RMA(dataT,within,{0.05,{'hsd'}});
print_for_manuscript(ra)
%%
n = 0;
[within,dvn,xlabels,awithinD] = make_within_table({'CN','TN'},[3,10]);
dataT = make_between_table({avar},dvn);
ra = RMA(dataT,within,{0.05,{''}});
ra.ranova
print_for_manuscript(ra)
%%
redF = [2]; redV = {2};
[dataTR,withinR] = reduce_within_between(dataT,within,redF,redV);
raR = RMA(dataTR,withinR,{0.05,{'hsd'}});
print_for_manuscript(raR)
%%
% visualization
mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size; dcolors = mData.dcolors;
tcolors = repmat(mData.dcolors(1),1,2);
% figure(300);clf; ha = gca;
ff = makeFigureRowsCols(2020,[10 4 1.25 1],'RowsCols',[1 1],'spaceRowsCols',[0.07 0],'rightUpShifts',[0.3 0.35],'widthHeightAdjustment',[-550 -380]);
MY = 25; ysp = 3.25; mY = 0; ystf = 3.12; ysigf = 0.05;titletxt = ''; ylabeltxt = {'PDF'}; % for all cells (vals) MY = 80
[hbs,xdata,mVar,semVar,combs,p,h] = view_results_rmanova(ff.h_axes(1,1),ra,{'AP','hsd',0.05},[1 2],tcolors,[mY MY ysp ystf ysigf],mData);
make_bars_hollow(hbs(2))
format_axes(gca);
set(gca,'xcolor','k','ycolor','k','xlim',xlim,'ylim',ylim,...
    'XTick',xdata,'XTickLabel',{'Air-On','Air-Off'});xtickangle(20);
ylabel({'Avg. Speed','(cm/s)'});
set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,2,{'Pooled'},{[0 0]});
% ht = set_axes_top_text_no_line(ff.hf,gca,sprintf('C1 - AOn'),[0.051 0.0 0 0]); 
save_pdf(ff.hf,mData.pdf_folder,sprintf('bar_graphs.pdf'),600);

%%
% visualization
mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size; dcolors = mData.dcolors;
tcolors = repmat(mData.dcolors(1),1,2);
% figure(300);clf; ha = gca;
ff = makeFigureRowsCols(2020,[10 4 1.25 1],'RowsCols',[1 1],'spaceRowsCols',[0.07 0],'rightUpShifts',[0.3 0.35],'widthHeightAdjustment',[-550 -380]);
MY = 3; ysp = 0.5; mY = -1; ystf = 0.52; ysigf = 0.05;titletxt = ''; ylabeltxt = {'PDF'}; % for all cells (vals) MY = 80
[hbs,xdata,mVar,semVar,combs,p,h] = view_results_rmanova(ff.h_axes(1,1),ra,{'AP','hsd',0.05},[1 2],tcolors,[mY MY ysp ystf ysigf],mData);
make_bars_hollow(hbs(2))
format_axes(gca);
set(gca,'xcolor','k','ycolor','k','xlim',xlim,'ylim',ylim,...
    'XTick',xdata,'XTickLabel',{'Air-On','Air-Off'});xtickangle(20);
ylabel({'Ske. Speed'});
set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,2,{'Pooled'},{[0 0]});
% ht = set_axes_top_text_no_line(ff.hf,gca,sprintf('C1 - AOn'),[0.051 0.0 0 0]); 
save_pdf(ff.hf,mData.pdf_folder,sprintf('bar_graphs.pdf'),600);




%%
% visualization
mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size; dcolors = mData.dcolors;
tcolors = repmat(mData.dcolors(1:10),1,2);
% figure(300);clf; ha = gca;
ff = makeFigureRowsCols(2020,[10 4 3.5 1],'RowsCols',[1 1],'spaceRowsCols',[0.07 0],'rightUpShifts',[0.07 0.35],'widthHeightAdjustment',[-100 -380]);
MY = 5; ysp = 3.25; mY = -1.5; ystf = 3.12; ysigf = 0.05;titletxt = ''; ylabeltxt = {'PDF'}; % for all cells (vals) MY = 80
[hbs,xdata,mVar,semVar,combs,p,h] = view_results_rmanova(ff.h_axes(1,1),ra,{'AP','hsd',0.05},[1 2],tcolors,[mY MY ysp ystf ysigf],mData);
format_axes(gca);
set(gca,'xcolor','k','ycolor','k','xlim',xlim,'ylim',ylim,...
    'XTick',xdata,'XTickLabel',{'Time-B','Dist-B'});xtickangle(20);
ylabel('Time (s)')
% set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,10,{'Time-Bin','Dist-Bin'},{[0 0]});
% ht = set_axes_top_text_no_line(ff.hf,gca,sprintf('C1 - AOn'),[0.051 0.0 0 0]); 
save_pdf(ff.hf,mData.pdf_folder,sprintf('bar_graphs.pdf'),600);