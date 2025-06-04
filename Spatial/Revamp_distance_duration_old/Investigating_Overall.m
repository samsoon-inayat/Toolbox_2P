function Investigating_Overall
%%
%Listing all the variables that I have related to the animal I have may be 
% . Time, Distance, Speed, and neural
%firing rates.
% experimental configurations include C3, C4, and C5 as well as presence or
% absence of stimuli include the air on and off phases and also the light
% stimulus in C4. Now first just to quantify motion, I can look at the time
% distance and speed and how they are changing with respect to air and
% light stimuli. This would be to quantify motion. The repeated measures
% ANOVA at the animal level will provide me with the information whether
% there was a certain effect of categorical variables which may be
% consistent across animals if I find significant differences.

udataT = evalin('base','udataT'); %Time Bins
udataD = evalin('base','udataD'); %Dist Bins
ei = evalin('base','ei');
configurations = {'C3','C4','C5'};
air_phases = {'ON','OFF'};
bin_types = {'time_bin','dist_bin'};
variable_combs = {'time_dist','time_speed','dist_speed','FR_time','FR_dist','FR_speed'};
variable_combs = {'time_dist','time_speed','dist_speed'};
azMI_vals = []; azPC_vals = [];
ow = [0 0] ;
for an = 1:5
    zMI_vals = []; zPC_vals = [];
    for bn = 1:length(bin_types)
        if bn == 1
            data_an = udataT{an};
        else
            data_an = udataD{an};
        end
        for cn = 1:length(configurations)
            for ap = 1:length(air_phases)
                out = get_continuous_variables(data_an,air_phases{ap},configurations{cn}); % concatenate the variables, time, distance, speed, and neural firing rates for trials in air phase (on or off) and configuration (C3, C4, or C5)
                % find the pairwise metrics including mutual information and
                % correlation
                for vn = 1:length(variable_combs)
                    tvar = variable_combs{vn};
                    spos = strfind(tvar,'_'); var1 = tvar(1:(spos-1)); var2 = tvar((spos+1):end);
                    cmdTxt = sprintf('var1v = out.%s;',var1); eval(cmdTxt); cmdTxt = sprintf('var2v = out.%s;',var2); eval(cmdTxt);
                    params = {'no_of_bins_for_MI',10,'no_of_shuffles_for_norm',500,'animal_info',ei{an},'overwrite_processing',ow,'air_phase',air_phases{ap},...
                        'configuration',configurations{cn},'variables',variable_combs{vn},'bin_type',bin_types{bn},'trial_type','concatenated'};
                    met = myMetrics(var1v,var2v,params);
                    MIs = nanmean(met.MI,1); PCs = nanmean(met.PC,1);
                    zMI_vals = [zMI_vals MIs(1,1)]; zPC_vals = [zPC_vals PCs(1,1)];
                end
            end
        end
    end
    azMI_vals(an,:) = zMI_vals; azPC_vals(an,:) = zPC_vals;
end

n = 0;
% (Intercept):Bin_Type:Air_Phase:Corr_Type [F(2,8) = 28.96, p < 0.001, η2 = .14] <--

%%
[within,dvn,xlabels,awithinD] = make_within_table({'BT','CN','AP','PT'},[2,3,2,3]);
data_matrix = azPC_vals; 
dataT = make_between_table({data_matrix},dvn);
ra = RMA(dataT,within,{0.05,{''}});
ra.ranova
print_for_manuscript(ra)
%% pooling data across configurations
% Initialize a matrix to store pooled data (5 animals × 12 conditions)
[dataP,withinP,withinPD] = pool_within_between(ra,'CN');
[within,dvn,xlabels,awithinD] = make_within_table({'BT','AP','PT'},[2,2,3]);
dataTP = make_between_table({dataP},dvn);
raP = RMA(dataTP,within,{0.05,{''}});
print_for_manuscript(raP)
%%
redF = [3]; redV = {3};
[dataTR1,withinR1] = reduce_within_between(dataTP,within,redF,redV);
raR = RMA(dataTR1,withinR1,{0.05/3,{'hsd'}});
print_for_manuscript(raR)
%%
% visualization
mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size; dcolors = mData.dcolors;
tcolors = repmat(mData.dcolors,1,10);
% figure(300);clf; ha = gca;
ff = makeFigureRowsCols(2020,[10 4 1 1],'RowsCols',[1 1],'spaceRowsCols',[0.07 0],'rightUpShifts',[0.27 0.25],'widthHeightAdjustment',[-360 -360]);
MY = 0.7; ysp = 0.15; mY = -0.5; ystf = 0.12; ysigf = 0.05;titletxt = ''; ylabeltxt = {'PDF'}; % for all cells (vals) MY = 80
[hbs,xdata,mVar,semVar,combs,p,h] = view_results_rmanova(ff.h_axes(1,1),raR,{'AP','hsd',(0.05/3)},[1 2],tcolors,[mY MY ysp ystf ysigf],mData);
format_axes(gca);
set(gca,'xcolor','k','ycolor','k','xlim',xlim,'ylim',ylim,...
    'XTick',xdata,'XTickLabel',{'AOn','AOff'});xtickangle(20);
ylabel('PCC')
box off; ht = set_axes_top_text_no_line(ff.hf,gca,sprintf('Dist-Speed'),[0 0.10 0.2 0]); 
% set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,2,{'Time-Bin','Dist-Bin'},{[0 0]});
save_pdf(ff.hf,mData.pdf_folder,sprintf('bar_graphs.pdf'),600);
%%
redF = [2]; redV = {1};
[dataTR2,withinR2] = reduce_within_between(dataTR1,withinR1,redF,redV);
raRR = RMA(dataTR2,withinR2,{(0.05/3)/2,{'hsd'}});
print_for_manuscript(raRR)

% visualization
mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size; dcolors = mData.dcolors;
tcolors = repmat(mData.dcolors,1,10);
% figure(300);clf; ha = gca;
ff = makeFigureRowsCols(2020,[10 4 1 1],'RowsCols',[1 1],'spaceRowsCols',[0.07 0],'rightUpShifts',[0.25 0.25],'widthHeightAdjustment',[-360 -360]);
MY = 1.5; ysp = 0.35; mY = 0; ystf = 0.2; ysigf = 0.05;titletxt = ''; ylabeltxt = {'PDF'}; % for all cells (vals) MY = 80
[hbs,xdata,mVar,semVar,combs,p,h] = view_results_rmanova(ff.h_axes(1,1),raRR,{'BT','hsd',(0.05/3)/2},[1 2],tcolors,[mY MY ysp ystf ysigf],mData);
format_axes(gca);
set(gca,'xcolor','k','ycolor','k','xlim',xlim,'ylim',ylim,...
    'XTick',xdata,'XTickLabel',{'Time-Bin','Dist-Bin'});xtickangle(20);
ylabel('PCC')
box off; ht = set_axes_top_text_no_line(ff.hf,gca,sprintf('Time-Distance (Air-ON)'),[0 0.10 0.2 0]); 
% set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,2,{'Time-Bin','Dist-Bin'},{[0 0]});
save_pdf(ff.hf,mData.pdf_folder,sprintf('bar_graphs.pdf'),600);
%% MI Values
[within,dvn,xlabels,awithinD] = make_within_table({'BT','CN','AP','PT'},[2,3,2,3]);
data_matrix = azMI_vals;
dataT = make_between_table({data_matrix},dvn);
ra = RMA(dataT,within,{0.05,{''}});
ra.ranova
print_for_manuscript(ra)

%%
redF = [4]; redV = {3};
[dataTR,withinR] = reduce_within_between(dataT,within,redF,redV);
raR = RMA(dataTR,withinR,{0.05,{'hsd'}});
print_for_manuscript(raR)

%%
% visualization
mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size; dcolors = mData.dcolors;
tcolors = repmat(mData.dcolors(3:5),1,2);
% figure(300);clf; ha = gca;
ff = makeFigureRowsCols(2020,[12 4 1.5 1],'RowsCols',[1 1],'spaceRowsCols',[0.07 0],'rightUpShifts',[0.15 0.3],'widthHeightAdjustment',[-230 -360]);
MY = 5; ysp = 0.35; mY = 0; ystf = 0.532; ysigf = 0.05;titletxt = ''; ylabeltxt = {'PDF'}; % for all cells (vals) MY = 80
[hbs,xdata,mVar,semVar,combs,p,h] = view_results_rmanova(ff.h_axes(1,1),ra,{'BT:CN','hsd',0.05},[1 2],tcolors,[mY MY ysp ystf ysigf],mData);
format_axes(gca);
set(gca,'xcolor','k','ycolor','k','xlim',xlim,'ylim',ylim,...
    'XTick',xdata,'XTickLabel',{'C3','C4','C5'});xtickangle(20);
ylabel('MIV')
box off; ht = set_axes_top_text_no_line(ff.hf,gca,sprintf('Time-Distance (Air-ON)'),[0 0.10 0.2 0]); 
set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,3,{['Time-Bin'],'Dist-Bin'},{[0 0]});
save_pdf(ff.hf,mData.pdf_folder,sprintf('bar_graphs.pdf'),600);

%%
% visualization
mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size; dcolors = mData.dcolors;
tcolors = repmat(mData.dcolors,1,10);
% figure(300);clf; ha = gca;
ff = makeFigureRowsCols(2020,[10 4 1.5 1],'RowsCols',[1 1],'spaceRowsCols',[0.07 0],'rightUpShifts',[0.25 0.3],'widthHeightAdjustment',[-360 -360]);
MY = 5; ysp = 0.35; mY = 0; ystf = 0.2; ysigf = 0.05;titletxt = ''; ylabeltxt = {'PDF'}; % for all cells (vals) MY = 80
[hbs,xdata,mVar,semVar,combs,p,h] = view_results_rmanova(ff.h_axes(1,1),ra,{'BT:AP','hsd',0.05},[1 2],tcolors,[mY MY ysp ystf ysigf],mData);
format_axes(gca);
set(gca,'xcolor','k','ycolor','k','xlim',xlim,'ylim',ylim,...
    'XTick',xdata,'XTickLabel',{'AOn','AOff'});xtickangle(20);
ylabel('MIV')
% box off; ht = set_axes_top_text_no_line(ff.hf,gca,sprintf('Time-Distance (Air-ON)'),[0 0.10 0.2 0]); 
set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,2,{['Time-' ...
    'Bin'],'Dist-Bin'},{[0 0]});
save_pdf(ff.hf,mData.pdf_folder,sprintf('bar_graphs.pdf'),600);

%% plot rasters for certain events
an = 1;
data_an = udataT{an};
ts = data_an.ts;

% air and C3 will give us configuration 3 and the air events but you have
% to detect the rising edges and falling edges. From the rising edges and
% falling edges you can easily find the air on phases by simple finding the
% variable of interest from a certain rising edge to a falling edge. Then
% the air off phase would be from a falling edge to a rising edge but for
% the last one, you have to extend based on the time;
csel = data_an.air & data_an.C3;
figure(100);clf;plot(ts,csel);pause(0.05);
redges = find_rising_edge(csel,0.5,0.05);
fedges = find_falling_edge(csel,-0.5,0.05);

if length(redges) == 10 && length(fedges) == 10
    disp('good')
end

sevent = redges; eevent = fedges; % for air on
sevent = fedges; eevent = [redges(2:end) ceil(fedges(end)+(mean(redges(2:end)-fedges(1:(length(fedges)-1)))*(ts(2)-ts(1))))]; % for air off
%% Overall Speed Tuning
% from the whole traces during the no-brake case, I will find out overall
% speed tuning by binning the speed and finding the average firing rate
% corresponding to it

% 1. find the nb and represent as a binary variable
for an = 1:5
    data_an = udataT{an};
    redges = find_rising_edge(data_an.bnb,0.5,0.05);
    fedges = find_falling_edge(data_an.bnb,-0.5,0.05);
    data_an.nb = logical(zeros(size(data_an.bnb)));
    data_an.nb((fedges(1)-1):(redges(2)+1))=1;
    udataT{an} = data_an;
end
%% run the speed tuning code to get firing_rate vs the speed (bins)
tic
for an = 1:5
    data_an = udataT{an};
    firing_rate = data_an.firing_rate;
    speed = data_an.speed;
    no_brake = data_an.nb;
    firing_rate = firing_rate(:,no_brake);
    speed = speed(no_brake);
    speed_tuning(an) = find_cell_speed_tuning(firing_rate,speed);
end
toc
%% calculate the mutual information (z-scored) for firing rates vs speed bins
nofmibins = 10; nshuffle = 500;
for an = 1:5
    FR_vs_speed = speed_tuning(an).FR_vs_speed;
    clear frs_MIz
    parfor ii = 1:size(FR_vs_speed,1)
        [output ~] = info_metrics_S(FR_vs_speed(ii,:), [], nofmibins, [], nshuffle); frs_MIz(ii) = output;
    end
    speed_tuning(an).MIz = frs_MIz;
end
%% for each animal make a list of speed cells from good gaussian fitting
for an = 1:5
    speed_struct = speed_tuning(an); temp_struct = speed_tuning(an).MIz;
    fieldname = 'ShannonMI_p'; spval = arrayfun(@(s) s.(fieldname), temp_struct);
    fieldname = 'corr_p'; cpval = arrayfun(@(s) s.(fieldname), temp_struct);
    bcs = speed_struct.bin_centers;
    fitg = speed_struct.fits.gauss;
    [rs,MFR,centers,PWs] = get_gauss_fit_parameters(fitg.coeffsrs,bcs(2)-bcs(1));
    inds = ~(centers < 1 | centers > 39 | rs < 0.25 | PWs < 10);% | PWs > 20 | PWs < 10;
    speed_tuning(an).good_gauss = inds;
    fraction_speed_cells(an) = sum(inds)/length(inds);
    fraction_speed_cells_cp(an) = sum(cpval < 0.05)/length(inds);
end
disp('Done');
beep;
%% this code is to view the individual graphs for cells
% we can view for selected cells defined in the variable inds. Firing rate
% vs speed bins is plotted
an = 3;
out = speed_tuning(an);
speed_struct = speed_tuning(an); temp_struct = speed_tuning(an).MIz;
fieldname = 'ShannonMI_p'; spval = arrayfun(@(s) s.(fieldname), temp_struct);
fieldname = 'corr_p'; cpval = arrayfun(@(s) s.(fieldname), temp_struct);
while 1
    d.bcs = out.bin_centers;
    d.FR = out.FR_vs_speed;
    fitg = out.fits.gauss; fits = out.fits.sigmoid; fitl = out.fits.linear;
    d.fFRl = fitl.fitted; d.fFRs = fits.fitted; d.fFRg = fitg.fitted;
    d.cl = fitl.coeffsrs(:,3); d.cs = fits.coeffsrs(:,3); d.cg = fitg.coeffsrs(:,3);
    [rs,MFR,centers,PWs] = get_gauss_fit_parameters(fitg.coeffsrs,d.bcs(2)-d.bcs(1));
    inds = ~(centers < 1 | centers > 39 | rs < 0.25 | PWs < 10);% | PWs > 20 | PWs < 10;
    inds = ~(rs < 0.25 | PWs < 10);% | PWs > 20 | PWs < 10;
    inds = cpval < 0.05 & inds;
    % inds = inds3 | inds4 | inds5;
    % inds = logical(resp_speed{an,4});
%     t_resp = cell_list_op(resp,[],'or');
%     inds = ~inds;
%     inds = inds & ~t_resp{an}';
    100*sum(inds)/length(inds)
    d.FR = d.FR(inds,:); d.fFRl = d.fFRl(inds,:); d.fFRs = d.fFRs(inds,:); d.fFRg = d.fFRg(inds,:);
    d.cl = d.cl(inds); d.cs = d.cs(inds); d.cg = d.cg(inds);
    generic_multi_plot(1000,[3,4,size(d.FR,1)],'plotSpeedTuning',d)
    break;
end


%% run the acceleration tuning code to get firing_rate vs the speed (bins)
tic
for an = 1:5
    data_an = udataT{an};
    firing_rate = data_an.firing_rate;
    speed = data_an.speed;
    ts = data_an.ts;
    ac = [0 diff(speed)./diff(ts)];
    no_brake = data_an.nb;
    firing_rate = firing_rate(:,no_brake);
    ac = ac(no_brake);
    accel_tuning(an) = find_cell_speed_tuning(firing_rate,ac);
end
toc
beep
%% calculate the mutual information (z-scored) for firing rates vs accel bins
tic
nofmibins = 10; nshuffle = 500;
for an = 1:5
    FR_vs_accel = accel_tuning(an).FR_vs_speed;
    clear frs_MIz
    parfor ii = 1:size(FR_vs_accel,1)
        [output ~] = info_metrics_S(FR_vs_accel(ii,:), [], nofmibins, [], nshuffle); frs_MIz(ii) = output;
    end
    accel_tuning(an).MIz = frs_MIz;
end
toc
beep;
%% for each animal make a list of speed cells from good gaussian fitting

for an = 1:5
    speed_struct = accel_tuning(an); temp_struct = accel_tuning(an).MIz;
    fieldname = 'ShannonMI_p'; spval = arrayfun(@(s) s.(fieldname), temp_struct);
    fieldname = 'corr_p'; cpval = arrayfun(@(s) s.(fieldname), temp_struct);
    bcs = speed_struct.bin_centers;
    fitg = speed_struct.fits.gauss;
    [rs,MFR,centers,PWs] = get_gauss_fit_parameters(fitg.coeffsrs,bcs(2)-bcs(1));
    inds = ~(centers < 1 | centers > 39 | rs < 0.25 | PWs < 10);% | PWs > 20 | PWs < 10;
    accel_tuning(an).good_gauss = inds;
    fraction_accel_cells(an) = sum(inds)/length(inds);
    fraction_accel_cells_cp(an) = sum(cpval < 0.05)/length(inds);
end
disp('Done');beep;
%% this code is to view the individual graphs for cells
% we can view for selected cells defined in the variable inds. Firing rate
% vs speed bins is plotted
an = 2;
out = accel_tuning(an);
speed_struct = accel_tuning(an); temp_struct = accel_tuning(an).MIz;
fieldname = 'ShannonMI_p'; spval = arrayfun(@(s) s.(fieldname), temp_struct);
fieldname = 'corr_p'; cpval = arrayfun(@(s) s.(fieldname), temp_struct);
while 1
    d.bcs = out.bin_centers;
    d.FR = out.FR_vs_speed;
    fitg = out.fits.gauss; fits = out.fits.sigmoid; fitl = out.fits.linear;
    d.fFRl = fitl.fitted; d.fFRs = fits.fitted; d.fFRg = fitg.fitted;
    d.cl = fitl.coeffsrs(:,3); d.cs = fits.coeffsrs(:,3); d.cg = fitg.coeffsrs(:,3);
    [rs,MFR,centers,PWs] = get_gauss_fit_parameters(fitg.coeffsrs,d.bcs(2)-d.bcs(1));
    inds = ~(centers < 1 | centers > 39 | rs < 0.25 | PWs < 10);% | PWs > 20 | PWs < 10;
    inds = ~(rs < 0.25 | PWs < 10);% | PWs > 20 | PWs < 10;
    inds = cpval < 0.05;% & inds;
    % inds = inds3 | inds4 | inds5;
    % inds = logical(resp_speed{an,4});
%     t_resp = cell_list_op(resp,[],'or');
%     inds = ~inds;
%     inds = inds & ~t_resp{an}';
    100*sum(inds)/length(inds)
    d.FR = d.FR(inds,:); d.fFRl = d.fFRl(inds,:); d.fFRs = d.fFRs(inds,:); d.fFRg = d.fFRg(inds,:);
    d.cl = d.cl(inds); d.cs = d.cs(inds); d.cg = d.cg(inds);
    generic_multi_plot(1000,[3,4,size(d.FR,1)],'plotSpeedTuning',d)
    break;
end
%% Let us compare the mutual information of cells for time vs distance 
% tuning for speed vs non-speed cells
for an = 1:5
    sp_cells = speed_tuning(an).good_gauss;
    nsp_cells = ~sp_cells;
    fetch_zMI = cell2mat(props_C.zMI(an,1:6));
    sp_cells_zMI = nanmean(fetch_zMI(sp_cells,:));
    nsp_cells_zMI = nanmean(fetch_zMI(nsp_cells,:));
    zMIs_an(an,:) = [sp_cells_zMI nsp_cells_zMI];
end

[within,dvn,xlabels,awithinD] = make_within_table({'SP','DT','Conf'},[2,2,3]);
dataT = make_between_table({zMIs_an},dvn);
ra = RMA(dataT,within,{0.05,{''}});
ra.ranova
print_for_manuscript(ra)

tcolors = repmat(mData.colors,1,10);
figure(300);clf; ha = gca;
MY = 4; ysp = 0.23; mY = 0; ystf = 0.25; ysigf = 0.05;titletxt = ''; ylabeltxt = {'PDF'}; % for all cells (vals) MY = 80
view_results_rmanova(ha,ra,{'SP','hsd',0.05},[1 2],tcolors,[mY MY ysp ystf ysigf],mData);

%% GLM
an = 1;
data_an = udataT{an};
% tRs = Rs_C{an,1};
ts = data_an.ts;ds = data_an.ds;sp = data_an.speed;tm = data_an.animal_motion;
ac = diff(sp)./diff(ts); ac = [0 ac];
csel = data_an.air & data_an.C5;
figure(100);clf;plot(ts,csel);pause(0.05);
redges = find_rising_edge(csel,0.5,0.05);
fedges = find_falling_edge(csel,-0.5,0.05);
cell_num = 14;
cell_raster = []; lencr = [];
atts = []; atds = []; atsp = []; attm = []; aFR = []; atac = [];
otts = [];
for ii = 1:length(redges)
    trial = redges(ii):fedges(ii);
    tts = ts(trial)-ts(trial(1)); tds = ds(trial)-ds(trial(1)); tsp = sp(trial); ttm = tm(trial); tac = ac(trial)
    FR = data_an.firing_rate(:,trial);
    atts = [atts tts]; atds = [atds tds]; atsp = [atsp tsp]; attm = [attm ttm]; aFR = [aFR FR]; atac = [atac tac];
    otts = [otts ts(trial)];
end
% find cells which have a lot of zeros
figure(100);clf;plot(FR');

% Creating the model for firing rate as a function of speed, distance, time, and motion
X = [atsp', atds', atts', attm', atac'];    % Predictor matrix (continuous variables)
model_results = struct();
tic
for ii = 1:size(aFR,1)
    ii
    neuron_idx = ii;
    y = aFR(ii,:);  y = y';% Flatten firing rate matrix for analysis
    % Fit a GLM to the data
    mdl = fitglm(X, y, 'linear');   % This will return a model object with coefficients and other information
    model_results(neuron_idx).mdl = mdl;
    
    % Store specific pieces of information from mdl (e.g., coefficients)
    model_results(neuron_idx).coefficients = mdl.Coefficients.Estimate;
    model_results(neuron_idx).p_values = mdl.Coefficients.pValue;
    
    % Optionally, store more details, e.g., R-squared, residuals
    model_results(neuron_idx).rsq = mdl.Rsquared.Ordinary;
    model_results(neuron_idx).residuals = mdl.Residuals.Raw;
end
toc;
beep;
%%
% Initialize an array to store results
significant_neurons = struct();

% Loop through each neuron
for neuron_idx = 1:num_neurons
    mdl = model_results(neuron_idx).mdl;  % Access model for the current neuron
    
    % Get p-values for each predictor (Speed, Distance, Time, Acceleration, Tissue Motion)
    p_values = mdl.Coefficients.pValue(2:end);  % Skip intercept
    
    % Get R-squared for the model
    rsq = model_results(neuron_idx).rsq;
    
    % Check if any p-value is below the significance threshold (e.g., 0.05)
    significant_predictors = p_values < 0.05;  % True/False for significance
    significant_neurons(neuron_idx).p_values = p_values;
    significant_neurons(neuron_idx).significant_predictors = significant_predictors;
    significant_neurons(neuron_idx).rsq = rsq;
end

%%
firing_rate_matrix = aFR;
% Example: Plot the firing rate of neurons that are significantly tuned to speed
for neuron_idx = 1:num_neurons
    if significant_neurons(neuron_idx).significant_predictors(2)  % Checking if speed is significant
        % Plot firing rate vs speed
        figure(100);clf;
        scatter(atsp, firing_rate_matrix(neuron_idx, :));  % Assuming you have speed and firing_rate_matrix
        title(['Neuron ', num2str(neuron_idx), ' Firing Rate vs Speed']);
        xlabel('Speed (cm/s)');
        ylabel('Firing Rate');
        pause(0.1);
    end
end


%%
% Initialize a matrix to store significant predictors (1 = significant, 0 = not significant)
num_neurons = length(significant_neurons);  % Number of neurons
num_variables = 5;  % Number of predictor variables: speed, distance, time, acceleration, tissue motion
significant_matrix = zeros(num_neurons, num_variables);

% Loop through each neuron and assign significance based on p-values
for neuron_idx = 1:num_neurons
    % Access the significant predictors for the current neuron
    significant_matrix(neuron_idx, :) = significant_neurons(neuron_idx).significant_predictors;
end

% Display the significant predictor matrix
disp(significant_matrix);



ts = data_an.ts;ds = data_an.ds;sp = data_an.speed;tm = data_an.animal_motion; 
ac = diff(sp)./diff(ts); ac = [0 ac];
cmdTxt = sprintf('csel = data_an.air & data_an.%s;',conf); eval(cmdTxt);
redges = find_rising_edge(csel,0.5,0.05);
fedges = find_falling_edge(csel,-0.5,0.05); 
if length(fedges) == 11 && fedges(1) == 1
    fedges(1) = [];
end
atts = []; atds = []; atsp = []; attm = []; aFR = []; atac = [];otts = [];
% figure(100);clf;plot(ts,csel);pause(0.05);
% air = 0;
if strcmp(air,'ON')
    sevent = redges; eevent = fedges; 
    for ii = 1:length(redges)
        trial = sevent(ii):eevent(ii);
        tts = ts(trial)-ts(trial(1)); tds = ds(trial)-ds(trial(1)); tsp = sp(trial); ttm = tm(trial); tac = ac(trial);
        FR = data_an.firing_rate(:,trial);
        atts = [atts tts]; atds = [atds tds]; atsp = [atsp tsp]; attm = [attm ttm]; aFR = [aFR FR]; atac = [atac tac];
        otts = [otts ts(trial)];
    end
end
if strcmp(air,'OFF')
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
% X = [atts', atsp', atds'];%, atac', attm'];    % Predictor matrix (continuous variables)
out.time = atts';
out.dist = atds';
out.speed = atsp';
out.FR = aFR';