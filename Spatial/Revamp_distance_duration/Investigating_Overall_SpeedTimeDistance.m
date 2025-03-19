function Investigating_Overall_SpeedTimeDistance
mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size;
ei = evalin('base','ei');
udata = evalin('base','udata');
udata1 = evalin('base','udata1');
udataT = evalin('base','udataT');
udataD = evalin('base','udataD');

%% 
udata = evalin('base','udataT'); %Time Bins
td_corr_T = [];
for an = 1:5
    td_corr = [];
    data_an = udata{an};
    ctrl = 0;
    temp = get_time_distance_corr(data_an,1,3,ctrl); td_corr = [td_corr temp];    temp = get_time_distance_corr(data_an,1,4,ctrl); td_corr = [td_corr temp];    temp = get_time_distance_corr(data_an,1,5,ctrl); td_corr = [td_corr temp];
    ctrl = 0;
    temp = get_time_distance_corr(data_an,0,3,ctrl); td_corr = [td_corr temp];    temp = get_time_distance_corr(data_an,0,4,ctrl); td_corr = [td_corr temp];    temp = get_time_distance_corr(data_an,0,5,ctrl); td_corr = [td_corr temp];
    td_corr_T = [td_corr_T;td_corr];
end

udata = evalin('base','udataD'); %Dist Bins
td_corr_D = [];
for an = 1:5
    td_corr = [];
    data_an = udata{an};
    ctrl = 0;
    temp = get_time_distance_corr(data_an,1,3,ctrl); td_corr = [td_corr temp];    temp = get_time_distance_corr(data_an,1,4,ctrl); td_corr = [td_corr temp];    temp = get_time_distance_corr(data_an,1,5,ctrl); td_corr = [td_corr temp];
    ctrl = 0;
    temp = get_time_distance_corr(data_an,0,3,ctrl); td_corr = [td_corr temp];    temp = get_time_distance_corr(data_an,0,4,ctrl); td_corr = [td_corr temp];    temp = get_time_distance_corr(data_an,0,5,ctrl); td_corr = [td_corr temp];
    td_corr_D = [td_corr_D;td_corr];
end
%%
[within,dvn,xlabels,awithinD] = make_within_table({'Bin_Type','Air_Phase','Conf_Num','Corr_Type'},[2,2,3,3]);
data_matrix = [td_corr_T td_corr_D];
dataT = make_between_table({data_matrix},dvn);
ra = RMA(dataT,within,{0.05,{''}});
ra.ranova
print_for_manuscript(ra)

%% pooling data across configurations
% Initialize a matrix to store pooled data (5 animals × 12 conditions)
[dataP,withinP,withinPD] = pool_within_between(ra,'Conf_Num');
[within,dvn,xlabels,awithinD] = make_within_table({'Bin_Type','Air_Phase','Corr_Type'},[2,2,3]);
dataTP = make_between_table({dataP},dvn);
raP = RMA(dataT,within,{0.05,{''}});
raP.ranova
print_for_manuscript(raP)
%%
redF = [3]; redV = {3};
[dataTR,withinR] = reduce_within_between(dataTP,withinP,redF,redV);
raR = RMA(dataTR,withinR,{0.05/3,{'hsd'}});
raR.ranova
print_for_manuscript(raR)
%%
tcolors = repmat(mData.colors,1,10);
figure(300);clf; ha = gca;
MY = 4; ysp = 0.23; mY = 0; ystf = 0.25; ysigf = 0.05;titletxt = ''; ylabeltxt = {'PDF'}; % for all cells (vals) MY = 80
view_results_rmanova(ha,ra,{'TD:Air','hsd',0.05},[1 2],tcolors,[mY MY ysp ystf ysigf],mData);


n = 0;



function [td_corr] = get_time_distance_corr(data_an,air,conf,ctrl)
ts = data_an.ts;ds = data_an.ds;sp = data_an.speed;tm = data_an.animal_motion; 

ac = diff(sp)./diff(ts); ac = [0 ac];
cmdTxt = sprintf('csel = data_an.air & data_an.C%d;',conf); eval(cmdTxt);
redges = find_rising_edge(csel,0.5,0.05);
fedges = find_falling_edge(csel,-0.5,0.05); 
if length(fedges) == 11 && fedges(1) == 1
    fedges(1) = [];
end
atts = []; atds = []; atsp = []; attm = []; aFR = []; atac = [];otts = [];
% figure(100);clf;plot(ts,csel);pause(0.05);
% air = 0;
if air 
    sevent = redges; eevent = fedges; 
    for ii = 1:length(redges)
        trial = sevent(ii):eevent(ii);
        tts = ts(trial)-ts(trial(1)); tds = ds(trial)-ds(trial(1)); tsp = sp(trial); ttm = tm(trial); tac = ac(trial);
        FR = data_an.firing_rate(:,trial);
        atts = [atts tts]; atds = [atds tds]; atsp = [atsp tsp]; attm = [attm ttm]; aFR = [aFR FR]; atac = [atac tac];
        otts = [otts ts(trial)];
    end
else
    sevent = fedges; eevent = [redges(2:end) ceil(fedges(end)+mean(redges(2:end)-fedges(1:(length(fedges)-1))))]; % for air off
    for ii = 1:length(redges)
        trial = sevent(ii):eevent(ii);
        tts = ts(trial)-ts(trial(1)); tds = ds(trial)-ds(trial(1)); tsp = sp(trial); ttm = tm(trial); tac = ac(trial);
        FR = data_an.firing_rate(:,trial);
        atts = [atts tts]; atds = [atds tds]; atsp = [atsp tsp]; attm = [attm ttm]; aFR = [aFR FR]; atac = [atac tac];
        otts = [otts ts(trial)];
    end
end

% Creating the model for firing rate as a function of speed, distance, time, and motion
X = [atts', atsp', atds'];%, atac', attm'];    % Predictor matrix (continuous variables)
corr_matrix = corr(X);
% disp('Correlation Matrix:');
% disp(corr_matrix);
td_corr = [corr_matrix(1,3),corr_matrix(1,2),corr_matrix(2,3)];
%% make a plot
if ctrl == 1 || ctrl == -1
    %
    mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size; dcolors = mData.dcolors;
    ff = makeFigureRowsCols(2020,[10 4 1.7 1.8],'RowsCols',[2 1],'spaceRowsCols',[0.07 0],'rightUpShifts',[0.15 0.15],'widthHeightAdjustment',[-360 -160]);
    timecol = dcolors{1}; distcol = dcolors{2}, speecol = dcolors{3};
    main_ax = ff.h_axes(1,1); axes(main_ax)
    [ax, h1, h2] = plotyy(1:length(atts),atts,1:length(atts),atds);set(h1,'color',timecol);set(h2,'color',distcol);
    ylabel(ax(1),'Time (s)'); ylabel(ax(2),'Distance (cm)');
    % xlabel('Time Bins');
    gca = ax(1);format_axes(gca);set(gca,'xcolor',colors{5},'ycolor',timecol,'xlim',[1 length(atts)/2],'ylim',[0 20],'XTickLabel',[]);
    gca = ax(2);format_axes(gca);set(gca,'xcolor','k','ycolor',distcol,'xlim',[1 length(atts)/2],'ylim',[0 150],'XTickLabel',[]);
    if ctrl == 1
        box off; ht = set_axes_top_text_no_line(ff.hf,main_ax,sprintf('C3 - Air-On (Concatenated Trials)'),[-0.1 0.080 0.2 0]); set(ht,'color',dcolors{7})
    end
    if ctrl == -1
        box off; ht = set_axes_top_text_no_line(ff.hf,main_ax,sprintf('C3 - Air-Off (Concatenated Trials)'),[-0.1 0.080 0.2 0]); set(ht,'color',dcolors{7})
    end
    ht = set_axes_top_text_no_line(ff.hf,ax(1),sprintf('r = %.3f',corr_matrix(1,3)),[0.25 0 0 0]); set(ht,'color',colors{7})

    main_ax = ff.h_axes(2,1); axes(main_ax)
    [ax, h1, h2] = plotyy(1:length(atts),atts,1:length(atsp),atsp);set(h1,'color',timecol);set(h2,'color',speecol);
    % ylabel(ax(1),'Distance (cm)'); ylabel(ax(2),'Speed (cm/s');
    ylabel(ax(1),'Time (s)'); ylabel(ax(2),'Speed (cm/s)');
    xlabel('Time Bins');
    gca = ax(1);format_axes(gca);set(gca,'xcolor',colors{5},'ycolor',timecol,'xlim',[1 length(atts)/2],'ylim',[0 20]);
    gca = ax(2);format_axes(gca);set(gca,'xcolor','k','ycolor',speecol,'xlim',[1 length(atts)/2],'ylim',[0 50]);
    box off;ht = set_axes_top_text_no_line(ff.hf,ax(1),sprintf('r = %.3f',corr_matrix(2,3)),[0.25 0 0 0]); set(ht,'color',colors{7})
    save_pdf(ff.hf,mData.pdf_folder,sprintf('time_dist_sp__bin.pdf'),600);
end
n = 0;
%% make a plot
if ctrl == 2 || ctrl == -2
    %
    mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size; dcolors = mData.dcolors;
    ff = makeFigureRowsCols(2020,[10 4 1.7 1.8],'RowsCols',[2 1],'spaceRowsCols',[0.07 0],'rightUpShifts',[0.15 0.15],'widthHeightAdjustment',[-360 -160]);
    xlimfac = 1.5;
    timecol = dcolors{1}; distcol = dcolors{2}, speecol = dcolors{3};
    main_ax = ff.h_axes(1,1); axes(main_ax)
    [ax, h1, h2] = plotyy(1:length(atts),atts,1:length(atts),atds);set(h1,'color',timecol);set(h2,'color',distcol);
    ylabel(ax(1),'Time (s)'); ylabel(ax(2),'Distance (cm)');
    % xlabel('Time Bins');
    gca = ax(1);format_axes(gca);set(gca,'xcolor',colors{5},'ycolor',timecol,'xlim',[1 length(atts)/xlimfac],'ylim',[0 20],'XTickLabel',[]);
    gca = ax(2);format_axes(gca);set(gca,'xcolor','k','ycolor',distcol,'xlim',[1 length(atts)/xlimfac],'ylim',[0 150],'XTickLabel',[]);
    if ctrl == 2
        box off; ht = set_axes_top_text_no_line(ff.hf,main_ax,sprintf('C3 - Air-On (Concatenated Trials)'),[-0.1 0.080 0.2 0]); set(ht,'color',dcolors{7})
    end
    if ctrl == -2
        box off; ht = set_axes_top_text_no_line(ff.hf,main_ax,sprintf('C3 - Air-Off (Concatenated Trials)'),[-0.1 0.080 0.2 0]); set(ht,'color',dcolors{7})
    end
    ht = set_axes_top_text_no_line(ff.hf,ax(1),sprintf('r = %.3f',corr_matrix(1,3)),[0.25 0 0 0]); set(ht,'color',colors{7})

    main_ax = ff.h_axes(2,1); axes(main_ax)
    [ax, h1, h2] = plotyy(1:length(atts),atts,1:length(atsp),atsp);set(h1,'color',timecol);set(h2,'color',speecol);
    % ylabel(ax(1),'Distance (cm)'); ylabel(ax(2),'Speed (cm/s');
    ylabel(ax(1),'Time (s)'); ylabel(ax(2),'Speed (cm/s)');
    xlabel('Distance Bins');
    gca = ax(1);format_axes(gca);set(gca,'xcolor',colors{5},'ycolor',timecol,'xlim',[1 length(atts)/xlimfac],'ylim',[0 20]);
    gca = ax(2);format_axes(gca);set(gca,'xcolor','k','ycolor',speecol,'xlim',[1 length(atts)/xlimfac],'ylim',[0 50]);
    box off;ht = set_axes_top_text_no_line(ff.hf,ax(1),sprintf('r = %.3f',corr_matrix(2,3)),[0.25 0 0 0]); set(ht,'color',colors{7})
    save_pdf(ff.hf,mData.pdf_folder,sprintf('time_dist_sp__bin.pdf'),600);
end
n = 0;

function find_pairwise_characteristics
%%
temp = 0;


% model_results = struct();
% tic
% for ii = 1:size(aFR,1)
%     ii
%     neuron_idx = ii;
%     y = aFR(ii,:);  y = y';% Flatten firing rate matrix for analysis
%     % Fit a GLM to the data
%     mdl = fitglm(X, y, 'linear');   % This will return a model object with coefficients and other information
%     model_results(neuron_idx).mdl = mdl;
% 
%     % Store specific pieces of information from mdl (e.g., coefficients)
%     model_results(neuron_idx).coefficients = mdl.Coefficients.Estimate;
%     model_results(neuron_idx).p_values = mdl.Coefficients.pValue;
% 
%     % Optionally, store more details, e.g., R-squared, residuals
%     model_results(neuron_idx).rsq = mdl.Rsquared.Ordinary;
%     model_results(neuron_idx).residuals = mdl.Residuals.Raw;
% end
% toc;
% beep;
% %%
% % Initialize an array to store results
% significant_neurons = struct();
% 
% % Loop through each neuron
% for neuron_idx = 1:num_neurons
%     mdl = model_results(neuron_idx).mdl;  % Access model for the current neuron
% 
%     % Get p-values for each predictor (Speed, Distance, Time, Acceleration, Tissue Motion)
%     p_values = mdl.Coefficients.pValue(2:end);  % Skip intercept
% 
%     % Get R-squared for the model
%     rsq = model_results(neuron_idx).rsq;
% 
%     % Check if any p-value is below the significance threshold (e.g., 0.05)
%     significant_predictors = p_values < 0.05;  % True/False for significance
%     significant_neurons(neuron_idx).p_values = p_values;
%     significant_neurons(neuron_idx).significant_predictors = significant_predictors;
%     significant_neurons(neuron_idx).rsq = rsq;
% end

%%
% firing_rate_matrix = aFR;
% % Example: Plot the firing rate of neurons that are significantly tuned to speed
% for neuron_idx = 1:num_neurons
%     if significant_neurons(neuron_idx).significant_predictors(2)  % Checking if speed is significant
%         % Plot firing rate vs speed
%         figure(100);clf;
%         scatter(atsp, firing_rate_matrix(neuron_idx, :));  % Assuming you have speed and firing_rate_matrix
%         title(['Neuron ', num2str(neuron_idx), ' Firing Rate vs Speed']);
%         xlabel('Speed (cm/s)');
%         ylabel('Firing Rate');
%         pause(0.1);
%     end
% end
% 
% 
% %%
% % Initialize a matrix to store significant predictors (1 = significant, 0 = not significant)
% num_neurons = length(significant_neurons);  % Number of neurons
% num_variables = 5;  % Number of predictor variables: speed, distance, time, acceleration, tissue motion
% significant_matrix = zeros(num_neurons, num_variables);
% 
% % Loop through each neuron and assign significance based on p-values
% for neuron_idx = 1:num_neurons
%     % Access the significant predictors for the current neuron
%     significant_matrix(neuron_idx, :) = significant_neurons(neuron_idx).significant_predictors;
% end
% 
% % Display the significant predictor matrix
% disp(significant_matrix);
