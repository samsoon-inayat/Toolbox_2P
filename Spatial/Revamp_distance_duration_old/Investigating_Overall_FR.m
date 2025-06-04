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
variable_combs = {'FR_time','FR_dist','FR_speed'};
aMIs = cell(5,36); aPCs = aMIs; aSigs = cell(5,36);
ow = [0 0] ;
for an = 1:5
    MI_vals = []; PC_vals = []; sigs = [];
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
                    met = myMetrics(var1v,var2v,params); MIs = met.MI; PCs = met.PC;
                    valsMI = {MIs}; valsPC = {PCs};
                    MI_vals = [MI_vals valsMI]; PC_vals = [PC_vals valsPC]; sigs = [sigs {out}];
                end
            end
        end
    end
    aMIs(an,:) = MI_vals; aPCs(an,:) = PC_vals; aSigs(an,:) = sigs;
end

n = 0;
% (Intercept):Bin_Type:Air_Phase:Corr_Type [F(2,8) = 28.96, p < 0.001, η2 = .14] <--
%%
% Above, I put a break in the code to get one animal data and then look at
% the tuning of different cells that are tuned to time, distance, or speed.
an = 3;
MIs = aMIs{an,3}; out = aSigs{an,3};
tuned_cells = MIs(:,3)<0.05 & ~isnan(MIs(:,1));
sum(tuned_cells)/length(tuned_cells)
% raster for a given cell
% Sample Data (Replace with your actual data)
% subs = [1; 1; 1; 2; 2; 2; 3; 3; 3; 3];  % Trial numbers
% vals = [5; 10; 15; 20; 25; 30; 35; 40; 45; 50]; % Firing rates
% time_bins = [1; 2; 3; 1; 2; 3; 1; 2; 3; 4]; % Corresponding time bins
cellnums = find(tuned_cells);
hf = figure(100);clf;
for ii = 1:length(cellnums)
    subs = out.trialnum;
    vals = out.FR(:,cellnums(ii));
    time_bins = out.binnum;
    
    % Determine number of trials and time bins
    num_trials = max(subs);
    num_time_bins = max(time_bins);
    
    % Create Raster Matrix: Trials (Rows) × Time Bins (Columns)
    raster_matrix = accumarray([subs, time_bins], vals, [num_trials, num_time_bins], @mean, NaN);
    
    % Create Raster Plot
    imagesc(raster_matrix); % Plot with imagesc
    colormap('hot'); % Set colormap for visualization
    colorbar; % Add colorbar for firing rate intensity
    xlabel('Time Bin'); 
    ylabel('Trial Number');
    title(sprintf('Raster Plot of Firing Rate across Trials and Time - %d',cellnums(ii)));
    
    % Improve Visualization
    set(gca, 'YDir', 'normal'); % Flip Y-axis to keep Trial 1 on top
    pause(0.051);
    k = get(gcf, 'CurrentCharacter'); % Get the last key pressed
    if ~isempty(k) && k == 27  % 27 is the ASCII code for Esc key
        disp('Esc key pressed. Exiting loop...');
        break;
    end
end
close(hf)

%%
[within,dvn,xlabels,awithinD] = make_within_table({'BT','CN','AP','PT'},[2,3,2,3]);
data_matrix = azMI_vals;
dataT = make_between_table({data_matrix},dvn);
ra = RMA(dataT,within,{0.05,{''}});
ra.ranova
print_for_manuscript(ra)

%%
redF = [2,3,4]; redV = {[3],[1],[1,2,3]};
redF = [1,3,4]; redV = {[1],[2],[1,2,3]};
[dataTR,withinR] = reduce_within_between(dataT,within,redF,redV);
raR = RMA(dataTR,withinR,{0.05,{'hsd'}});
print_for_manuscript(raR)

%
mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size; dcolors = mData.dcolors;
tcolors = repmat(mData.colors,1,10);
ff = makeFigureRowsCols(2020,[12 4 1.5 1],'RowsCols',[1 1],'spaceRowsCols',[0.07 0],'rightUpShifts',[0.15 0.3],'widthHeightAdjustment',[-230 -360]);
MY = 3; ysp = 0.1; mY = 0; ystf = 0.1; ysigf = 0.01;titletxt = ''; ylabeltxt = {'PDF'}; % for all cells (vals) MY = 80
[hbs,xdata,mVar,semVar,combs,p,h] = view_results_rmanova(ff.h_axes(1,1),raR,{'PT','hsd',0.05},[1 2],tcolors,[mY MY ysp ystf ysigf],mData);
format_axes(gca);
set(gca,'xcolor','k','ycolor','k','xlim',xlim,'ylim',ylim,...
    'XTick',xdata,'XTickLabel',{'Time','Dist','Speed'});xtickangle(20);
ylabel('MIV')
% box off; ht = set_axes_top_text_no_line(ff.hf,gca,sprintf('Time-Distance (Air-ON)'),[0 0.10 0.2 0]); 
set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,length(hbs)/2,{['Time-Bin'],'Dist-Bin'},{[0 0]});






