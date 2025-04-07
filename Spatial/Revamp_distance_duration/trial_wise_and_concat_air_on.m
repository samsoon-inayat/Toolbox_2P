function time_trial_distance_speed

udata = evalin('base','udata1');
ei = evalin('base','ei');
configurations = {'C3','C4','C5'};
air_phases = {'ON','OFF'};

bin_types = {'time_bin','dist_bin'};
variable_combs = {'time_dist','time_speed','dist_speed','FR_time','FR_dist','FR_speed'};
variable_combs = {'time_dist','time_speed','dist_speed'};
variable_combs = {'FR_time','FR_dist','FR_speed'};

time_bins = 0:0.3:1000; % set a large number of bins
dist_bins = 0:3:1000;
all_data = {};
for an = 1:5
    data_an = udata{an};
    field_names = fieldnames(data_an);
    for ii = 1:length(field_names)
        varname = field_names{ii};
        cmdTxt = sprintf('clear %s;',varname);eval(cmdTxt);
        cmdTxt = sprintf('%s = data_an.%s;',varname,varname);eval(cmdTxt);
    end
    firing_rate = data_an.firing_rate;
    frf_n = double(frf);
    num_frames = sum(frf==1);
    frf_n(frf==1) = 1:num_frames;
    % ds(ds<0) = NaN; var = fillmissing(ds, 'spline');
    % dds = diff(ds); spike_idx = find(abs(dds) > 50) + 1;
    % ds(spike_idx) = NaN; ds(spike_idx-1) = NaN; ds(spike_idx+1) = NaN; ds = fillmissing(ds, 'spline');

    for cn = 1:length(configurations)
        for ap = 1:2
            for tn = 1:10
                % get time to complete trial
                if ap == 1
                    cmdTxt = sprintf('air_trials = air_trials_on .* %s;',configurations{cn}); eval(cmdTxt);
                else
                    cmdTxt = sprintf('air_trials = air_trials_off .* %s;',configurations{cn}); eval(cmdTxt);
                end
                idx = find(air_trials == tn & frf);
                tts = ts(idx)-ts(idx(1));  
                tds = ds(idx)-ds(idx(1)); tds(tds<0) = 0;
                tsp = speed(idx);
                idx_fr = frf_n(idx);
                FR = firing_rate(:,idx_fr);
                for bt = 1:2
                    [an cn ap tn bt]
                    if bt == 1
                        bin_indices = discretize(tts,time_bins);
                    else
                        bin_indices = discretize(tds,dist_bins);
                    end
                    time_binned = accumarray(bin_indices',tts,[],@mean); time_binned = time_binned - time_binned(1);
                    dist_binned = accumarray(bin_indices',tds,[],@mean); dist_binned = dist_binned - dist_binned(1);
                    speed_binned = accumarray(bin_indices',tsp,[],@mean); 

                    FR_binned = [];
                    parfor neuron_idx = 1:size(firing_rate,1)
                        % Calculate the binned firing rate for each neuron using the frame-based bin indices
                        FR_binned(neuron_idx, :) = accumarray(bin_indices', FR(neuron_idx, :),[], @mean, NaN)';
                    end
                    FR_binned = FR_binned';
                    td.time_binned = time_binned; td.dist_binned = dist_binned; td.speed_binned = speed_binned; td.FR_binned = FR_binned;
                    all_data{an,cn,ap,tn,bt} = td;
                    n = 0;
                end
            end
        end
    end
end
n = 0;
%% concat trial-wise
var_names = {'time','dist','speed','FR'};
for an = 1:5
    for cn = 1:3
        for ap = 1:2
            for vn = 1:length(var_names)
                cmdTxt = sprintf('%s_t = [];',var_names{vn}); eval(cmdTxt);
                cmdTxt = sprintf('%s_d = [];',var_names{vn}); eval(cmdTxt);
            end
            for tn = 1:10
                for bt = 1:2
                    [an cn ap tn bt]
                    td = all_data{an,cn,ap,tn,bt};
                    if bt == 1
                        for vn = 1:length(var_names)
                            cmdTxt = sprintf('%s_t = [%s_t;td.%s_binned];',var_names{vn},var_names{vn},var_names{vn}); eval(cmdTxt);
                        end
                    else
                        for vn = 1:length(var_names)
                            cmdTxt = sprintf('%s_d = [%s_d;td.%s_binned];',var_names{vn},var_names{vn},var_names{vn}); eval(cmdTxt);
                        end
                    end
                end
            end
        end
    end
end

%%
nn = 1;
an = 1
anvarC = avarC{an};
anvar = []; 
for cn = 1:length(configurations)
    for ap = 1:2
        for tn = 1:10
            for bt = 1:2
                for vn = 1:length(variable_combs)
                    [an cn ap tn bt vn]
                    vidx = 3;
                    thisvar_vn1 = anvarC{cn,ap,tn,bt,1}; thisvar_vn2 = anvarC{cn,ap,tn,bt,2}; thisvar_vn3 = anvarC{cn,ap,tn,bt,3};
                    thisvar = [thisvar_vn1(:,vidx) thisvar_vn2(:,vidx) thisvar_vn3(:,vidx)];
                    anvar = [anvar thisvar(nn,1)];% outD.trial_metrics(:,idx)'];
                end
            end
        end
    end
end
avar = anvar;
%%
n = 0;
clc
[within,dvn,xlabels,awithinD] = make_within_table({'CN','AP','TN','BT','PT'},[3,2,10,2,3]);
dataT = make_between_table({avar},dvn);
ra = RMA(dataT,within,{0.05,{''}});
print_for_manuscript(ra)
%%
n = 0;
[within,dvn,xlabels,awithinD] = make_within_table({'CN','TN','BT'},[3,10,2]);
dataT = make_between_table({avar},dvn);
ra = RMA(dataT,within,{0.05,{''}});
print_for_manuscript(ra)
%%
clc
raR = RMA_bonferroni(ra,5);
%%
[xdata,mVar,semVar,combs,p,h] = get_vals_RMA(mData,raR{3},{'AP:BT','hsd',0.05},[1 2]);

%%
temp_data = table2array(raR{1}.rm.BetweenDesign);
descriptiveStatistics(temp_data(:));
%%
clc
raRR = RMA_bonferroni(raR{3},4);
%%
temp_data = table2array(raRR{1}.rm.BetweenDesign);
descriptiveStatistics(temp_data(:));
%%
n = 0;
[within,dvn,xlabels,awithinD] = make_within_table({'CN','TN'},[3,10]);
dataT = make_between_table({avar},dvn);
ra = RMA(dataT,within,{0.05,{''}});
ra.ranova
print_for_manuscript(ra)
%%
redF = [4,5]; redV = {1,1};
[dataTR,withinR] = reduce_within_between(dataT,within,redF,redV);
raR = RMA(dataTR,withinR,{0.05,{'hsd'}});
print_for_manuscript(raR)
%%
% visualization
mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size; dcolors = mData.dcolors;
tcolors = repmat(mData.dcolors(1:10),1,3);
% figure(300);clf; ha = gca;
ff = makeFigureRowsCols(2020,[4 4 6.25 1],'RowsCols',[1 1],'spaceRowsCols',[0.07 0],'rightUpShifts',[0.3 0.35],'widthHeightAdjustment',[-550 -380]);
MY = 0.5; ysp = 1.25; mY = 0; ystf = 1.12; ysigf = 0.05;titletxt = ''; ylabeltxt = {'PDF'}; % for all cells (vals) MY = 80
[hbs,xdata,mVar,semVar,combs,p,h] = view_results_rmanova(ff.h_axes(1,1),ra,{'CN:TN','hsd',0.05},[1 2],tcolors,[mY MY ysp ystf ysigf],mData);
% make_bars_hollow(hbs(2))
format_axes(gca);
xlbl = {'T01','T02','T03','T04','T05','T06','T07','T08','T09','T10'}
set(gca,'xcolor','k','ycolor','k','xlim',xlim,'ylim',ylim,...
    'XTick',xdata,'XTickLabel',xlbl);xtickangle(20);
ylabel({'Mutual','Information'});
% set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,2,{'Pooled'},{[0 0]});
ht = set_axes_top_text_no_line(ff.hf,gca,sprintf('C1 - AOn'),[0.051 0.0 0 0]); 
save_pdf(ff.hf,mData.pdf_folder,sprintf('bar_graphs.pdf'),600);

%%
% visualization
mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size; dcolors = mData.dcolors;
tcolors = repmat(mData.dcolors(1),1,2);
% figure(300);clf; ha = gca;
ff = makeFigureRowsCols(2020,[10 4 1.25 1],'RowsCols',[1 1],'spaceRowsCols',[0.07 0],'rightUpShifts',[0.3 0.35],'widthHeightAdjustment',[-550 -380]);
MY = 3; ysp = 0.5; mY = -1; ystf = 0.52; ysigf = 0.05;titletxt = ''; ylabeltxt = {'PDF'}; % for all cells (vals) MY = 80
[hbs,xdata,mVar,semVar,combs,p,h] = view_results_rmanova(ff.h_axes(1,1),raR{2},{'AP','hsd',0.05},[1 2],tcolors,[mY MY ysp ystf ysigf],mData);
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
MY = 1; ysp = 0.25; mY = -1; ystf = 0.12; ysigf = 0.05;titletxt = ''; ylabeltxt = {'PDF'}; % for all cells (vals) MY = 80
[hbs,xdata,mVar,semVar,combs,p,h] = view_results_rmanova(ff.h_axes(1,1),raR{3},{'AP:BT','hsd',0.05},[1 2],tcolors,[mY MY ysp ystf ysigf],mData);
format_axes(gca);
set(gca,'xcolor','k','ycolor','k','xlim',xlim,'ylim',ylim,...
    'XTick',xdata,'XTickLabel',{'Time-B','Dist-B'});xtickangle(20);
ylabel('MI')
% set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,10,{'Time-Bin','Dist-Bin'},{[0 0]});
% ht = set_axes_top_text_no_line(ff.hf,gca,sprintf('C1 - AOn'),[0.051 0.0 0 0]); 
save_pdf(ff.hf,mData.pdf_folder,sprintf('bar_graphs.pdf'),600);