function dist_dur
%% analysing the resp. RF, zMI, and RS for all cells
ntrials = 50;
C2_AOn_T = Ab_t_T; C2_AOff_T = Ab_i_T;  C7_AOn_T = Abs_t_T; C7_AOff_T = Abs_i_T;
C3_AOn_T = Ar_t_T; C3_AOff_T = Ar_i_T;   
C4_AOn_T = ArL_t_T; C4_AOff_T = ArL_i_T;
C5_AOn_T = Ars_t_T; C5_AOff_T = Ars_i_T;

C3_AOn_D = Ar_t_D; C3_AOff_D = Ar_i_D; 
C4_AOn_D = ArL_t_D; C4_AOff_D = ArL_i_D;
C5_AOn_D = Ars_t_D; C5_AOff_D = Ars_i_D;

si_old = [Ab_t_T Ab_i_T Ar_t_T Ar_i_T ArL_t_T ArL_i_T Ars_t_T Ars_i_T Ar_t_D Ar_i_D ArL_t_D ArL_i_D Ars_t_D Ars_i_D Abs_t_T Abs_i_T];
si1 = [C2_AOn_T C2_AOff_T C3_AOn_T C3_AOff_T C4_AOn_T C4_AOff_T C5_AOn_T C5_AOff_T C3_AOn_D C3_AOff_D C4_AOn_D C4_AOff_D C5_AOn_D C5_AOff_D C7_AOn_T C7_AOff_T];
%     si = [C1_i_T C2_i_T C3_i_T C4_i_T];
Rs_C = o.Rs(:,si1);mRs_C = o.mR(:,si1);
props_C = get_props_Rs(Rs_C,ntrials);
pop_var_name = {'all','vals','valsT','Nvals','good_zMI','Ngood_zMI'};
pop_var_name = {'all'};
% pop_var_name = {'good_FR'};
sel_pop_C = cell_list_op(props_C,pop_var_name);

si = {'C2_AOn_T','C2_AOff_T','C3_AOn_T','C3_AOff_T','C4_AOn_T','C4_AOff_T','C5_AOn_T','C5_AOff_T','C3_AOn_D','C3_AOff_D',...
    'C4_AOn_D','C4_AOff_D','C5_AOn_D','C5_AOff_D','C7_AOn_T','C7_AOff_T'};

zMI_cell = props_C.zMI;
num_animals = 5;
%% extract zMI values into a big matrix for C2 and C7 and their phases
% Define the relevant indices from the 'si' sequence for C2 and C7 time tuning
% Assuming the sequence 'si' holds the conditions in the exact order you provided:
C2_AOn_T_idx = find(strcmp(si, 'C2_AOn_T'));  % Index for C2 Air-On Time Tuning
C2_AOff_T_idx = find(strcmp(si, 'C2_AOff_T'));  % Index for C2 Air-Off Time Tuning
C7_AOn_T_idx = find(strcmp(si, 'C7_AOn_T'));  % Index for C7 Air-On Time Tuning
C7_AOff_T_idx = find(strcmp(si, 'C7_AOff_T'));  % Index for C7 Air-Off Time Tuning

% Initialize a cell array to hold zMI matrices for each animal
zMI_brake_matrices = cell(num_animals, 1);

% Loop through each animal to create the zMI matrix
for ani_idx = 1:num_animals
    % Extract the zMI vectors for each condition (M x 1 vectors)
    zMI_C2_AOn_T = zMI_cell{ani_idx, C2_AOn_T_idx};  % C2 Air-On Time Tuning
    zMI_C2_AOff_T = zMI_cell{ani_idx, C2_AOff_T_idx};  % C2 Air-Off Time Tuning
    zMI_C7_AOn_T = zMI_cell{ani_idx, C7_AOn_T_idx};  % C7 Air-On Time Tuning
    zMI_C7_AOff_T = zMI_cell{ani_idx, C7_AOff_T_idx};  % C7 Air-Off Time Tuning
    
    % Get the number of neurons (M) for this animal, assume zMI vectors have the same length
    num_neurons = size(zMI_C2_AOn_T, 1); 
    
    % Initialize the zMI matrix for this animal (M neurons x 4 conditions)
    zMI_matrix_animal = NaN(num_neurons, 4);
    
    % Fill the zMI matrix with the extracted zMI values
    zMI_matrix_animal(:, 1) = zMI_C2_AOn_T;   % C2 Air-On Time Tuning
    zMI_matrix_animal(:, 2) = zMI_C2_AOff_T;  % C2 Air-Off Time Tuning
    zMI_matrix_animal(:, 3) = zMI_C7_AOn_T;   % C7 Air-On Time Tuning
    zMI_matrix_animal(:, 4) = zMI_C7_AOff_T;  % C7 Air-Off Time Tuning
    
    % Store the matrix for this animal
    zMI_brake_matrices{ani_idx} = zMI_matrix_animal;
end

% Now, zMI_brake_matrices holds the zMI matrix for each animal with 4 columns:
% [C2_AOn_T, C2_AOff_T, C7_AOn_T, C7_AOff_T]


%% extract zMI values into a big matrix for C3, C4, and C5 and their phases (time tuning)

% Define the relevant indices from the 'si' sequence for C3, C4, C5 time tuning
C3_AOn_T_idx = find(strcmp(si, 'C3_AOn_T'));  % Index for C3 Air-On Time Tuning
C3_AOff_T_idx = find(strcmp(si, 'C3_AOff_T'));  % Index for C3 Air-Off Time Tuning
C4_AOn_T_idx = find(strcmp(si, 'C4_AOn_T'));  % Index for C4 Air-On Time Tuning
C4_AOff_T_idx = find(strcmp(si, 'C4_AOff_T'));  % Index for C4 Air-Off Time Tuning
C5_AOn_T_idx = find(strcmp(si, 'C5_AOn_T'));  % Index for C5 Air-On Time Tuning
C5_AOff_T_idx = find(strcmp(si, 'C5_AOff_T'));  % Index for C5 Air-Off Time Tuning

% Initialize a cell array to hold zMI matrices for each animal
zMI_no_brake_matrices = cell(num_animals, 1);

% Loop through each animal to create the zMI matrix
for ani_idx = 1:num_animals
    % Extract the zMI vectors for each condition (M x 1 vectors)
    zMI_C3_AOn_T = zMI_cell{ani_idx, C3_AOn_T_idx};  % C3 Air-On Time Tuning
    zMI_C3_AOff_T = zMI_cell{ani_idx, C3_AOff_T_idx};  % C3 Air-Off Time Tuning
    zMI_C4_AOn_T = zMI_cell{ani_idx, C4_AOn_T_idx};  % C4 Air-On Time Tuning
    zMI_C4_AOff_T = zMI_cell{ani_idx, C4_AOff_T_idx};  % C4 Air-Off Time Tuning
    zMI_C5_AOn_T = zMI_cell{ani_idx, C5_AOn_T_idx};  % C5 Air-On Time Tuning
    zMI_C5_AOff_T = zMI_cell{ani_idx, C5_AOff_T_idx};  % C5 Air-Off Time Tuning
    
    % Get the number of neurons (M) for this animal, assume zMI vectors have the same length
    num_neurons = size(zMI_C3_AOn_T, 1); 
    
    % Initialize the zMI matrix for this animal (M neurons x 6 conditions)
    zMI_matrix_animal = NaN(num_neurons, 6);
    
    % Fill the zMI matrix with the extracted zMI values
    zMI_matrix_animal(:, 1) = zMI_C3_AOn_T;   % C3 Air-On Time Tuning
    zMI_matrix_animal(:, 2) = zMI_C3_AOff_T;  % C3 Air-Off Time Tuning
    zMI_matrix_animal(:, 3) = zMI_C4_AOn_T;   % C4 Air-On Time Tuning
    zMI_matrix_animal(:, 4) = zMI_C4_AOff_T;  % C4 Air-Off Time Tuning
    zMI_matrix_animal(:, 5) = zMI_C5_AOn_T;   % C5 Air-On Time Tuning
    zMI_matrix_animal(:, 6) = zMI_C5_AOff_T;  % C5 Air-Off Time Tuning
    
    % Store the matrix for this animal
    zMI_no_brake_matrices{ani_idx} = zMI_matrix_animal;
end

% Now, zMI_no_brake_matrices holds the zMI matrix for each animal with 6 columns:
% [C3_AOn_T, C3_AOff_T, C4_AOn_T, C4_AOff_T, C5_AOn_T, C5_AOff_T]


%% extract zMI values into a big matrix for C3, C4, and C5 and their phases (distance tuning)
% Define the relevant indices from the 'si' sequence for C3, C4, C5 distance tuning
C3_AOn_D_idx = find(strcmp(si, 'C3_AOn_D'));  % Index for C3 Air-On Distance Tuning
C3_AOff_D_idx = find(strcmp(si, 'C3_AOff_D'));  % Index for C3 Air-Off Distance Tuning
C4_AOn_D_idx = find(strcmp(si, 'C4_AOn_D'));  % Index for C4 Air-On Distance Tuning
C4_AOff_D_idx = find(strcmp(si, 'C4_AOff_D'));  % Index for C4 Air-Off Distance Tuning
C5_AOn_D_idx = find(strcmp(si, 'C5_AOn_D'));  % Index for C5 Air-On Distance Tuning
C5_AOff_D_idx = find(strcmp(si, 'C5_AOff_D'));  % Index for C5 Air-Off Distance Tuning

% Initialize a cell array to hold zMI matrices for each animal
zMI_no_brake_dist_matrices = cell(num_animals, 1);

% Loop through each animal to create the zMI matrix
for ani_idx = 1:num_animals
    % Extract the zMI vectors for each condition (M x 1 vectors)
    zMI_C3_AOn_D = zMI_cell{ani_idx, C3_AOn_D_idx};  % C3 Air-On Distance Tuning
    zMI_C3_AOff_D = zMI_cell{ani_idx, C3_AOff_D_idx};  % C3 Air-Off Distance Tuning
    zMI_C4_AOn_D = zMI_cell{ani_idx, C4_AOn_D_idx};  % C4 Air-On Distance Tuning
    zMI_C4_AOff_D = zMI_cell{ani_idx, C4_AOff_D_idx};  % C4 Air-Off Distance Tuning
    zMI_C5_AOn_D = zMI_cell{ani_idx, C5_AOn_D_idx};  % C5 Air-On Distance Tuning
    zMI_C5_AOff_D = zMI_cell{ani_idx, C5_AOff_D_idx};  % C5 Air-Off Distance Tuning
    
    % Get the number of neurons (M) for this animal, assume zMI vectors have the same length
    num_neurons = size(zMI_C3_AOn_D, 1); 
    
    % Initialize the zMI matrix for this animal (M neurons x 6 conditions)
    zMI_matrix_animal = NaN(num_neurons, 6);
    
    % Fill the zMI matrix with the extracted zMI values
    zMI_matrix_animal(:, 1) = zMI_C3_AOn_D;   % C3 Air-On Distance Tuning
    zMI_matrix_animal(:, 2) = zMI_C3_AOff_D;  % C3 Air-Off Distance Tuning
    zMI_matrix_animal(:, 3) = zMI_C4_AOn_D;   % C4 Air-On Distance Tuning
    zMI_matrix_animal(:, 4) = zMI_C4_AOff_D;  % C4 Air-Off Distance Tuning
    zMI_matrix_animal(:, 5) = zMI_C5_AOn_D;   % C5 Air-On Distance Tuning
    zMI_matrix_animal(:, 6) = zMI_C5_AOff_D;  % C5 Air-Off Distance Tuning
    
    % Store the matrix for this animal
    zMI_no_brake_dist_matrices{ani_idx} = zMI_matrix_animal;
end

% Now, zMI_no_brake_dist_matrices holds the zMI matrix for each animal with 6 columns:
% [C3_AOn_D, C3_AOff_D, C4_AOn_D, C4_AOff_D, C5_AOn_D, C5_AOff_D]


%% mean for no brake conditions to perform ANOVA

% Assuming zMI_no_brake_matrices contains time tuning values
% And zMI_no_brake_dist_matrices contains distance tuning values

% Initialize matrices to store mean zMI values for each animal
mean_zMI_time = NaN(num_animals, 6);  % 6 columns for C3, C4, C5 air-on and air-off phases (time)
mean_zMI_dist = NaN(num_animals, 6);  % 6 columns for C3, C4, C5 air-on and air-off phases (distance)

% Loop over animals to calculate mean zMI for each configuration and phase
for ani_idx = 1:num_animals
    % Time tuning mean values over neurons
    mean_zMI_time(ani_idx, :) = mean(zMI_no_brake_matrices{ani_idx}, 1, 'omitnan');  % Average over neurons for each column (condition)
    
    % Distance tuning mean values over neurons
    mean_zMI_dist(ani_idx, :) = mean(zMI_no_brake_dist_matrices{ani_idx}, 1, 'omitnan');  % Average over neurons for each column (condition)
end

% Now, 'mean_zMI_time' and 'mean_zMI_dist' contain the average zMI values for each animal
% Rows = animals, Columns = [C3_AOn, C3_AOff, C4_AOn, C4_AOff, C5_AOn, C5_AOff]


[within,dvn,xlabels] = make_within_table({'E','C','P'},[2,3,2]);
dataT = make_between_table({[mean_zMI_time mean_zMI_dist]},dvn);
raNB = RMA(dataT,within,{0.05,{'hsd'}});
%     ra.ranova
print_for_manuscript(raNB);

%% mean for brake conditions to perform ANOVA
% Assuming zMI_brake_matrices contains time tuning values for C2 and C7

% Initialize a matrix to store mean zMI values for each animal
mean_zMI_brake_time = NaN(num_animals, 4);  % 4 columns for C2_AOn_T, C2_AOff_T, C7_AOn_T, C7_AOff_T (time tuning)

% Loop over animals to calculate mean zMI for each configuration and phase
for ani_idx = 1:num_animals
    % Time tuning mean values over neurons for C2 and C7
    mean_zMI_brake_time(ani_idx, :) = mean(zMI_brake_matrices{ani_idx}, 1, 'omitnan');  % Average over neurons for each column (condition)
end

% Now, 'mean_zMI_brake_time' contains the average zMI values for each animal
% Rows = animals, Columns = [C2_AOn_T, C2_AOff_T, C7_AOn_T, C7_AOff_T]

[within,dvn,xlabels] = make_within_table({'C','P'},[2,2]);
dataT = make_between_table({mean_zMI_brake_time},dvn);
raB = RMA(dataT,within,{0.05,{'hsd'}});
%     ra.ranova
print_for_manuscript(raB);

%% comparing across brake and no-brake condition

% For C2 and C7 (Brake condition) - From 'mean_zMI_brake_time'
mean_zMI_C2 = mean_zMI_brake_time(:, 1:2);  % C2_AOn_T, C2_AOff_T
mean_zMI_C7 = mean_zMI_brake_time(:, 3:4);  % C7_AOn_T, C7_AOff_T

% For C3 and C5 (No-Brake condition) - From 'mean_zMI_time'
mean_zMI_C3 = mean_zMI_time(:, 1:2);  % C3_AOn_T, C3_AOff_T
mean_zMI_C5 = mean_zMI_time(:, 5:6);  % C5_AOn_T, C5_AOff_T

[within,dvn,xlabels] = make_within_table({'B','P'},[2,2]);
dataT = make_between_table({[mean_zMI_C2 mean_zMI_C3]},dvn);
raC2C3 = RMA(dataT,within,{0.05,{'hsd'}});
%     ra.ranova
print_for_manuscript(raC2C3);


[within,dvn,xlabels] = make_within_table({'B','P'},[2,2]);
dataT = make_between_table({[mean_zMI_C7 mean_zMI_C5]},dvn);
raC7C5 = RMA(dataT,within,{0.05,{'hsd'}});
%     ra.ranova
print_for_manuscript(raC7C5);

%% mean over phases in each brake and no-brake conditions

% Brake Air-On (Average across C2_AOn_T and C7_AOn_T)
mean_zMI_brake_on = mean(mean_zMI_brake_time(:, [1, 3]), 2);  % Average C2_AOn_T and C7_AOn_T
% No-Brake Air-On (Average across C3_AOn_T, C4_AOn_T, C5_AOn_T)
mean_zMI_no_brake_on = mean(mean_zMI_time(:, [1, 3, 5]), 2);  % Average C3_AOn_T, C4_AOn_T, C5_AOn_T

% Brake Air-Off (Average across C2_AOff_T and C7_AOff_T)
mean_zMI_brake_off = mean(mean_zMI_brake_time(:, [1, 3]+1), 2);  % Average C2_AOn_T and C7_AOn_T
% No-Brake Air-Off (Average across C3_AOff_T, C4_AOff_T, C5_AOff_T)
mean_zMI_no_brake_off = mean(mean_zMI_time(:, [1, 3, 5]+1), 2);  % Average C3_AOn_T, C4_AOn_T, C5_AOn_T

[within,dvn,xlabels] = make_within_table({'B','P'},[2,2]);
dataT = make_between_table({[mean_zMI_brake_on mean_zMI_brake_off mean_zMI_no_brake_on mean_zMI_no_brake_off]},dvn);
raBnB = RMA(dataT,within,{0.05,{'hsd'}});
%     ra.ranova
print_for_manuscript(raBnB);


%% Figure zMI comparison for only Brake
print_for_manuscript(raB);
magfac = mData.magfac;
ff = makeFigureRowsCols(107,[3 5 2 1.25],'RowsCols',[1 1],'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.07 0.3],...
    'widthHeightAdjustment',[10 -450]);
MY = 0.5; ysp = 9; mY = 0; ystf = 7; ysigf = 1;titletxt = ''; ylabeltxt = {'PDF'}; % for all cells (vals) MY = 80
stp = 0.33*magfac; widths = [1.1 0.35 2.85 1]*magfac+0.061; gap = 0.09105*magfac;
adjust_axes(ff,[mY MY],stp,widths,gap,{''});
axes_title_shifts_line = [0.051 0.55 -0.1 0]; axes_title_shifts_text = [0.02 0.1 0 0]; xs_gaps = [1 1.5];

tcolors = repmat(mData.colors(1:2),1,6);
[hbs,xdata] = view_results_rmanova(ff.h_axes(1,1),raB,{'C:P','hsd',0.05},xs_gaps,tcolors,[mY MY ysp ystf ysigf],mData);
xtls = {'AOn','AOff'}; set(gca,'xtick',xdata,'xticklabels',xtls); xtickangle(30);
set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,2,{'C2','C7'},{[0.001 0.0051]});
ylabel('zMI');
axes_title(ff,{1},{'Brake Condition'},axes_title_shifts_line,axes_title_shifts_text);

save_pdf(ff.hf,mData.pdf_folder,sprintf('bar_graph.pdf'),600);
%% Figure zMI comparison for only No-Brake
print_for_manuscript(raNB);
magfac = mData.magfac;
ff = makeFigureRowsCols(107,[3 5 4.55 1.25],'RowsCols',[1 3],'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.07 0.4],...
    'widthHeightAdjustment',[10 -450]);
MY = 1.5; ysp = 0.31; mY = 0; ystf = 0.2; ysigf = 0.091;titletxt = ''; ylabeltxt = {'PDF'}; % for all cells (vals) MY = 80
stp = 0.33*magfac; widths = [3.1 0.35 0.35 1]*magfac+0.061; gap = 0.09105*magfac;
adjust_axes(ff,[mY MY],stp,widths,gap,{''});
axes_title_shifts_line = [0 0.55 0 0]; axes_title_shifts_text = [0.02 0.1 0 0]; xs_gaps = [1 1.5];

tcolors = repmat(mData.colors(1:2),1,6);
[hbs,xdata] = view_results_rmanova(ff.h_axes(1,1),raNB,{'E:C:P','hsd',0.05},xs_gaps,tcolors,[mY MY ysp ystf ysigf],mData);
xtls = {'AOn','AOff'}; set(gca,'xtick',xdata,'xticklabels',xtls); xtickangle(30);
set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,2,{'C3','C4','C5','C3','C4','C5'},{[0.001 0.0051]});
set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,6,{'Time','Dist'},{[-0.13 0.0051]});
ylabel('zMI');
% axes_title(ff,{1},{'Brake Condition'},axes_title_shifts_line,axes_title_shifts_text);


tcolors = repmat(mData.dcolors(3:4),1,6);
[hbs,xdata] = view_results_rmanova(ff.h_axes(1,2),raNB,{'E','hsd',0.05},xs_gaps,tcolors,[mY MY ysp ystf ysigf],mData);
xtls = {'Time','Dist'}; set(gca,'xtick',xdata,'xticklabels',xtls); xtickangle(30);
set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,2,{'Pooled'},{[0.001 0.0051]});

tcolors = repmat(mData.colors(1:end),1,6);
[hbs,xdata] = view_results_rmanova(ff.h_axes(1,3),raNB,{'P','hsd',0.05},xs_gaps,tcolors,[mY MY ysp ystf ysigf],mData);
xtls = {'AOn','AOff'}; set(gca,'xtick',xdata,'xticklabels',xtls); xtickangle(30);
set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,2,{'Pooled'},{[0.001 0.0051]});

save_pdf(ff.hf,mData.pdf_folder,sprintf('bar_graph.pdf'),600);


%% Figure zMI comparison for only Brake

magfac = mData.magfac;
ff = makeFigureRowsCols(107,[3 5 2.35 1.25],'RowsCols',[1 3],'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.07 0.4],...
    'widthHeightAdjustment',[10 -450]);
MY = 0.5; ysp = 0.1; mY = 0; ystf = 0.1; ysigf = 0.031;titletxt = ''; ylabeltxt = {'PDF'}; % for all cells (vals) MY = 80
stp = 0.33*magfac; widths = [0.85 0.35 0.35 1]*magfac+0.061; gap = 0.09105*magfac;
adjust_axes(ff,[mY MY],stp,widths,gap,{''});
axes_title_shifts_line = [0 0.55 0 0]; axes_title_shifts_text = [0.02 0.1 0 0]; xs_gaps = [1 1.5];

tcolors = repmat(mData.colors(1:2),1,2);
[hbs,xdata] = view_results_rmanova(ff.h_axes(1,1),raBnB,{'B:P','hsd',0.05,''},xs_gaps,tcolors,[mY MY ysp ystf ysigf],mData);
xtls = {'AOn','AOff'}; set(gca,'xtick',xdata,'xticklabels',xtls); xtickangle(30);
set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,2,{'Brake (Br)','No-Brake (NBr)'},{[0.001 0.0051]});ylabel('zMI');

tcolors = repmat(mData.dcolors(1:2),1,6);
[hbs,xdata] = view_results_rmanova(ff.h_axes(1,2),raBnB,{'B','hsd',0.05},xs_gaps,tcolors,[mY MY ysp ystf ysigf],mData);
xtls = {'Br','NBr'}; set(gca,'xtick',xdata,'xticklabels',xtls); xtickangle(30);
set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,2,{'Pooled'},{[0.001 0.0051]});

tcolors = repmat(mData.colors(1:end),1,6);
[hbs,xdata] = view_results_rmanova(ff.h_axes(1,3),raBnB,{'P','hsd',0.05},xs_gaps,tcolors,[mY MY ysp ystf ysigf],mData);
xtls = {'AOn','AOff'}; set(gca,'xtick',xdata,'xticklabels',xtls); xtickangle(30);
set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,2,{'Pooled'},{[0.001 0.0051]});

save_pdf(ff.hf,mData.pdf_folder,sprintf('bar_graph.pdf'),600);



% axes_title(ff,{1},{'Brake Condition'},axes_title_shifts_line,axes_title_shifts_text);


%% checking zMIs across phases and configurations for percent non-NaN values for Brake condition

% Assuming zMI_brake_matrices is a cell array storing data for each animal
num_animals = length(zMI_brake_matrices);  % Number of animals
num_brake_configs = 2;  % C2 and C7
num_phases = 2;  % Air-On and Air-Off

% Preallocate for the 5x4 matrix
percent_matrix = zeros(num_animals, num_brake_configs * num_phases);  % 5x4 matrix

% Loop through animals and configurations
for ani_idx = 1:num_animals
    col_idx = 1;  % Column index for the matrix
    for config_idx = 1:num_brake_configs
        % Extract zMI values for Brake condition (C2, C7)
        zMI_AOn = zMI_brake_matrices{ani_idx}(:, (config_idx - 1) * 2 + 1);  % Air-On
        zMI_AOff = zMI_brake_matrices{ani_idx}(:, (config_idx - 1) * 2 + 2);  % Air-Off
        
        % Calculate percentage of neurons with non-NaN zMI for Air-On and Air-Off
        percent_matrix(ani_idx, col_idx) = 100 * sum(~isnan(zMI_AOn)) / length(zMI_AOn);  % Air-On
        percent_matrix(ani_idx, col_idx + 1) = 100 * sum(~isnan(zMI_AOff)) / length(zMI_AOff);  % Air-Off
        
        % Increment column index for the next configuration
        col_idx = col_idx + 2;
    end
end

% Display the matrix
disp('5x4 matrix with percentage of non-NaN zMI values for C2 and C7, air-on and air-off:');
disp(percent_matrix);


% Create grouped bar plot for the 5x4 matrix
figure(1000); clf;
bar(percent_matrix);
xticks(1:num_animals);
% xtls({'Animal 1', 'Animal 2', 'Animal 3', 'Animal 4', 'Animal 5'});
xlabel('Animal');
ylabel('Percentage of Neurons with non-NaN zMI');
legend({'C2 Air-On', 'C2 Air-Off', 'C7 Air-On', 'C7 Air-Off'}, 'Location', 'NorthEastOutside');
title('Percentage of Neurons with non-NaN zMI for C2 and C7, Air-On and Air-Off Phases');


[within,dvn,xlabels] = make_within_table({'C','P'},[2,2]);
dataT = make_between_table({percent_matrix},dvn);
raB = RMA(dataT,within,{0.05,{'hsd'}});
%     ra.ranova
print_for_manuscript(raB);

%% Figure zMI comparison for only Brake
t_ra = raB;
print_for_manuscript(raB);
magfac = mData.magfac;
ff = makeFigureRowsCols(107,[3 5 2 1.25],'RowsCols',[1 1],'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.07 0.3],...
    'widthHeightAdjustment',[10 -450]);
MY = 115; ysp = 9; mY = 0; ystf = 7; ysigf = 1;titletxt = ''; ylabeltxt = {'PDF'}; % for all cells (vals) MY = 80
stp = 0.33*magfac; widths = [1.1 0.35 2.85 1]*magfac+0.061; gap = 0.09105*magfac;
adjust_axes(ff,[mY MY],stp,widths,gap,{''});
axes_title_shifts_line = [0.051 0.55 -0.1 0]; axes_title_shifts_text = [0.02 0.1 0 0]; xs_gaps = [1 1.5];

tcolors = repmat(mData.colors(1:2),1,6);
[hbs,xdata] = view_results_rmanova(ff.h_axes(1,1),t_ra,{'C:P','hsd',0.05},xs_gaps,tcolors,[mY MY ysp ystf ysigf],mData);
xtls = {'AOn','AOff'}; set(gca,'xtick',xdata,'xticklabels',xtls); xtickangle(30);
set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,2,{'C2','C7'},{[0.001 0.0051]});
ylabel('% non-NaN zMI Values');
axes_title(ff,{1},{'Brake Condition'},axes_title_shifts_line,axes_title_shifts_text);

save_pdf(ff.hf,mData.pdf_folder,sprintf('bar_graph.pdf'),600);


%% check non-NaN zMI values across configurations and phases during No-Brake condition
% Assuming zMI_no_brake_matrices stores time encoding data for each animal
num_animals = length(zMI_no_brake_matrices);  % Number of animals
num_no_brake_configs = 3;  % C3, C4, C5
num_phases = 2;  % Air-On and Air-Off

% Preallocate for the 5x6 matrix (3 configs x 2 phases)
percent_no_brake_time_matrix = zeros(num_animals, num_no_brake_configs * num_phases);  % 5x6 matrix

% Loop through animals and configurations for time encoding
for ani_idx = 1:num_animals
    col_idx = 1;  % Column index for the matrix
    for config_idx = 1:num_no_brake_configs
        % Extract zMI values for No-Brake condition (time encoding)
        zMI_AOn_time = zMI_no_brake_matrices{ani_idx}(:, (config_idx - 1) * 2 + 1);  % Air-On (Time)
        zMI_AOff_time = zMI_no_brake_matrices{ani_idx}(:, (config_idx - 1) * 2 + 2);  % Air-Off (Time)
        
        % Calculate percentage of neurons with non-NaN zMI for Air-On and Air-Off (Time)
        percent_no_brake_time_matrix(ani_idx, col_idx) = 100 * sum(~isnan(zMI_AOn_time)) / length(zMI_AOn_time);  % Air-On
        percent_no_brake_time_matrix(ani_idx, col_idx + 1) = 100 * sum(~isnan(zMI_AOff_time)) / length(zMI_AOff_time);  % Air-Off
        
        % Move to the next configuration
        col_idx = col_idx + 2;
    end
end



% Preallocate for the 5x6 matrix (3 configs x 2 phases)
percent_no_brake_dist_matrix = zeros(num_animals, num_no_brake_configs * num_phases);  % 5x6 matrix

% Loop through animals and configurations for distance encoding
for ani_idx = 1:num_animals
    col_idx = 1;  % Column index for the matrix
    for config_idx = 1:num_no_brake_configs
        % Extract zMI values for No-Brake condition (distance encoding)
        zMI_AOn_dist = zMI_no_brake_dist_matrices{ani_idx}(:, (config_idx - 1) * 2 + 1);  % Air-On (Distance)
        zMI_AOff_dist = zMI_no_brake_dist_matrices{ani_idx}(:, (config_idx - 1) * 2 + 2);  % Air-Off (Distance)
        
        % Calculate percentage of neurons with non-NaN zMI for Air-On and Air-Off (Distance)
        percent_no_brake_dist_matrix(ani_idx, col_idx) = 100 * sum(~isnan(zMI_AOn_dist)) / length(zMI_AOn_dist);  % Air-On
        percent_no_brake_dist_matrix(ani_idx, col_idx + 1) = 100 * sum(~isnan(zMI_AOff_dist)) / length(zMI_AOff_dist);  % Air-Off
        
        % Move to the next configuration
        col_idx = col_idx + 2;
    end
end


percent_no_brake_matrix = [percent_no_brake_time_matrix percent_no_brake_dist_matrix];

% Create grouped bar plot for the 5x12 matrix
figure(1002); clf;
bar(percent_no_brake_matrix);
xticks(1:num_animals);
xticklabels({'Animal 1', 'Animal 2', 'Animal 3', 'Animal 4', 'Animal 5'});
xlabel('Animal');
ylabel('Percentage of Neurons with non-NaN zMI');
legend({'C3 Time Air-On', 'C3 Time Air-Off', 'C4 Time Air-On', 'C4 Time Air-Off', 'C5 Time Air-On', 'C5 Time Air-Off', ...
    'C3 Distance Air-On', 'C3 Distance Air-Off', 'C4 Distance Air-On', 'C4 Distance Air-Off', 'C5 Distance Air-On', 'C5 Distance Air-Off'}, ...
    'Location', 'NorthEastOutside');
title('Percentage of Neurons with non-NaN zMI for C3, C4, and C5 (Time and Distance), Air-On and Air-Off Phases');


%
[within,dvn,xlabels] = make_within_table({'E','C','P'},[2,3,2]);
dataT = make_between_table({[percent_no_brake_matrix]},dvn);
raNB_P = RMA(dataT,within,{0.05,{'hsd'}});
%     ra.ranova
print_for_manuscript(raNB_P);

% %
% [within,dvn,xlabels] = make_within_table({'C','P'},[3,2]);
% dataT = make_between_table({[percent_no_brake_dist_matrix]},dvn);
% raNB_P = RMA(dataT,within,{0.05,{'hsd'}});
% %     ra.ranova
% print_for_manuscript(raNB_P);
% 
% %
% [within,dvn,xlabels] = make_within_table({'C','P'},[3,2]);
% dataT = make_between_table({[percent_no_brake_time_matrix]},dvn);
% raNB_P = RMA(dataT,within,{0.05,{'hsd'}});
% %     ra.ranova
% print_for_manuscript(raNB_P);

%% mean over phases in each brake and no-brake conditions

% Brake Air-On (Average across C2_AOn_T and C7_AOn_T)
mean_nzMI_brake_on = mean(percent_matrix(:, [1, 3]), 2);  % Average C2_AOn_T and C7_AOn_T
% No-Brake Air-On (Average across C3_AOn_T, C4_AOn_T, C5_AOn_T)
mean_nzMI_no_brake_on = mean(percent_no_brake_time_matrix(:, [1, 3, 5]), 2);  % Average C3_AOn_T, C4_AOn_T, C5_AOn_T

% Brake Air-Off (Average across C2_AOff_T and C7_AOff_T)
mean_nzMI_brake_off = mean(percent_matrix(:, [1, 3]+1), 2);  % Average C2_AOn_T and C7_AOn_T
% No-Brake Air-Off (Average across C3_AOff_T, C4_AOff_T, C5_AOff_T)
mean_nzMI_no_brake_off = mean(percent_no_brake_time_matrix(:, [1, 3, 5]+1), 2);  % Average C3_AOn_T, C4_AOn_T, C5_AOn_T

[within,dvn,xlabels] = make_within_table({'B','P'},[2,2]);
dataT = make_between_table({[mean_nzMI_brake_on mean_nzMI_brake_off mean_nzMI_no_brake_on mean_nzMI_no_brake_off]},dvn);
raBnB = RMA(dataT,within,{0.05,{'hsd'}});
%     ra.ranova
print_for_manuscript(raBnB);
%% Figure zMI comparison for only Brake

magfac = mData.magfac;
ff = makeFigureRowsCols(107,[3 5 2.35 1.25],'RowsCols',[1 1],'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.07 0.4],...
    'widthHeightAdjustment',[10 -450]);
MY = 160; ysp = 20; mY = 0; ystf = 40; ysigf = 5;titletxt = ''; ylabeltxt = {'PDF'}; % for all cells (vals) MY = 80
stp = 0.33*magfac; widths = [0.85 0.35 0.35 1]*magfac+0.061; gap = 0.09105*magfac;
adjust_axes(ff,[mY MY],stp,widths,gap,{''});
axes_title_shifts_line = [0 0.55 0 0]; axes_title_shifts_text = [0.02 0.1 0 0]; xs_gaps = [1 1.5];

tcolors = repmat(mData.colors(1:2),1,2);
[hbs,xdata] = view_results_rmanova(ff.h_axes(1,1),raBnB,{'B:P','hsd',0.05},xs_gaps,tcolors,[mY MY ysp ystf ysigf],mData);
xtls = {'AOn','AOff'}; set(gca,'xtick',xdata,'xticklabels',xtls); xtickangle(30);
set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,2,{'Brake (Br)','No-Brake (NBr)'},{[0.001 0.0051]});
ylabel('% non-NaN zMI Values');

% tcolors = repmat(mData.dcolors(1:2),1,6);
% [hbs,xdata] = view_results_rmanova(ff.h_axes(1,2),raBnB,{'B','hsd',0.05},xs_gaps,tcolors,[mY MY ysp ystf ysigf],mData);
% xtls = {'Br','NBr'}; set(gca,'xtick',xdata,'xticklabels',xtls); xtickangle(30);
% set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,2,{'Pooled'},{[0.001 0.0051]});
% 
% tcolors = repmat(mData.colors(1:end),1,6);
% [hbs,xdata] = view_results_rmanova(ff.h_axes(1,3),raBnB,{'P','hsd',0.05},xs_gaps,tcolors,[mY MY ysp ystf ysigf],mData);
% xtls = {'AOn','AOff'}; set(gca,'xtick',xdata,'xticklabels',xtls); xtickangle(30);
% set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,2,{'Pooled'},{[0.001 0.0051]});

save_pdf(ff.hf,mData.pdf_folder,sprintf('bar_graph.pdf'),600);



% axes_title(ff,{1},{'Brake Condition'},axes_title_shifts_line,axes_title_shifts_text);


%% Figure percent non-NaN zMI comparison for only No-Brake
t_ra = raNB_P;
print_for_manuscript(t_ra);
magfac = mData.magfac;
ff = makeFigureRowsCols(107,[3 5 4.55 1.25],'RowsCols',[1 3],'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.07 0.4],...
    'widthHeightAdjustment',[10 -450]);
MY = 125; ysp = 0.31; mY = 0; ystf = 0.2; ysigf = 0.091;titletxt = ''; ylabeltxt = {'PDF'}; % for all cells (vals) MY = 80
stp = 0.33*magfac; widths = [3.1 0.35 0.35 1]*magfac+0.061; gap = 0.09105*magfac;
adjust_axes(ff,[mY MY],stp,widths,gap,{''});
axes_title_shifts_line = [0 0.55 0 0]; axes_title_shifts_text = [0.02 0.1 0 0]; xs_gaps = [1 1.5];

tcolors = repmat(mData.colors(1:2),1,6);
[hbs,xdata] = view_results_rmanova(ff.h_axes(1,1),t_ra,{'E:C:P','hsd',0.05},xs_gaps,tcolors,[mY MY ysp ystf ysigf],mData);
xtls = {'AOn','AOff'}; set(gca,'xtick',xdata,'xticklabels',xtls); xtickangle(30);
set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,2,{'C3','C4','C5','C3','C4','C5'},{[0.001 0.0051]});
set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,6,{'Time','Dist'},{[-0.13 0.0051]});
ylabel('% non-NaN zMI Values');
% axes_title(ff,{1},{'Brake Condition'},axes_title_shifts_line,axes_title_shifts_text);


tcolors = repmat(mData.dcolors(5:7),1,6);
[hbs,xdata] = view_results_rmanova(ff.h_axes(1,2),t_ra,{'C','hsd',0.05},xs_gaps,tcolors,[mY MY ysp ystf ysigf],mData);
xtls = {'C3','C4','C5'}; set(gca,'xtick',xdata,'xticklabels',xtls); xtickangle(30);
set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,3,{'Pooled'},{[0.001 0.0051]});

tcolors = repmat(mData.colors(1:end),1,6);
[hbs,xdata] = view_results_rmanova(ff.h_axes(1,3),t_ra,{'P','hsd',0.05},xs_gaps,tcolors,[mY MY ysp ystf ysigf],mData);
xtls = {'AOn','AOff'}; set(gca,'xtick',xdata,'xticklabels',xtls); xtickangle(30);
set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,2,{'Pooled'},{[0.001 0.0051]});

save_pdf(ff.hf,mData.pdf_folder,sprintf('bar_graph.pdf'),600);

%%
% Preallocate for storing differences in zMI for neurons with non-NaN values
zMI_diff_time = cell(num_animals, num_no_brake_configs);  % For time encoding
zMI_diff_dist = cell(num_animals, num_no_brake_configs);  % For distance encoding

% Loop through animals and configurations
for ani_idx = 1:num_animals
    for config_idx = 1:num_no_brake_configs
        % Extract zMI values for Air-On and Air-Off (Time)
        zMI_AOn_time = zMI_no_brake_matrices{ani_idx}(:, (config_idx - 1) * 2 + 1);  % Air-On (Time)
        zMI_AOff_time = zMI_no_brake_matrices{ani_idx}(:, (config_idx - 1) * 2 + 2);  % Air-Off (Time)

        % Find neurons with non-NaN zMI values in both phases
        valid_neurons_time = ~isnan(zMI_AOn_time) & ~isnan(zMI_AOff_time);
        
        % Calculate the difference in zMI between phases (Air-On - Air-Off)
        zMI_diff_time{ani_idx, config_idx} = zMI_AOn_time(valid_neurons_time) - zMI_AOff_time(valid_neurons_time);
        
        % Extract zMI values for Air-On and Air-Off (Distance)
        zMI_AOn_dist = zMI_no_brake_dist_matrices{ani_idx}(:, (config_idx - 1) * 2 + 1);  % Air-On (Distance)
        zMI_AOff_dist = zMI_no_brake_dist_matrices{ani_idx}(:, (config_idx - 1) * 2 + 2);  % Air-Off (Distance)

        % Find neurons with non-NaN zMI values in both phases (Distance)
        valid_neurons_dist = ~isnan(zMI_AOn_dist) & ~isnan(zMI_AOff_dist);

        % Calculate the difference in zMI between phases (Air-On - Air-Off)
        zMI_diff_dist{ani_idx, config_idx} = zMI_AOn_dist(valid_neurons_dist) - zMI_AOff_dist(valid_neurons_dist);
    end
end
%%
% Plot histogram of zMI differences for a specific configuration (e.g., C3)
figure(1003); clf;
subplot(1,2,1);
histogram(zMI_diff_time{1, 1}, 'BinWidth', 0.1);  % Time encoding for Animal 1, C3
xlabel('zMI Difference (Air-On - Air-Off)');
ylabel('Number of Neurons');
title('Time Encoding: C3 zMI Differences Between Phases');

subplot(1,2,2);
histogram(zMI_diff_dist{1, 1}, 'BinWidth', 0.1);  % Distance encoding for Animal 1, C3
xlabel('zMI Difference (Air-On - Air-Off)');
ylabel('Number of Neurons');
title('Distance Encoding: C3 zMI Differences Between Phases');
% Aggregate zMI differences across all animals and configurations
zMI_diff_time_all = [];
zMI_diff_dist_all = [];

for ani_idx = 1:num_animals
    for config_idx = 1:num_no_brake_configs
        zMI_diff_time_all = [zMI_diff_time_all; zMI_diff_time{ani_idx, config_idx}];
        zMI_diff_dist_all = [zMI_diff_dist_all; zMI_diff_dist{ani_idx, config_idx}];
    end
end

% Plot aggregate histograms
figure(1004); clf;
subplot(1,2,1);
histogram(zMI_diff_time_all, 'BinWidth', 0.1);
xlabel('zMI Difference (Air-On - Air-Off)');
ylabel('Number of Neurons');
title('Time Encoding: Aggregate zMI Differences');

subplot(1,2,2);
histogram(zMI_diff_dist_all, 'BinWidth', 0.1);
xlabel('zMI Difference (Air-On - Air-Off)');
ylabel('Number of Neurons');
title('Distance Encoding: Aggregate zMI Differences');

%%
% Initialize a list to store the zMI differences for all animals and neurons
zMI_differences_C2 = [];
zMI_differences_C7 = [];

% Loop through each animal and calculate the zMI differences for C2 and C7
for ani_idx = 1:numel(zMI_brake_matrices)
    % Get zMI values for the current animal
    zMI_data = zMI_brake_matrices{ani_idx};
    zMI_data_nb = zMI_no_brake_matrices{ani_idx};
    
    % C2 zMI differences (Air-On - Air-Off)
    zMI_diff_C2 = zMI_data(:, 1) - zMI_data(:, 2);
    
    % C7 zMI differences (Air-On - Air-Off)
    zMI_diff_C7 = zMI_data(:, 3) - zMI_data(:, 4);
    
    % Append to the list of zMI differences
    zMI_differences_C2 = [zMI_differences_C2; zMI_diff_C2];
    zMI_differences_C7 = [zMI_differences_C7; zMI_diff_C7];
end

% Plot the distribution of zMI differences for C2
figure(1001); clf;
subplot(1, 2, 1);
histogram(zMI_differences_C2, 'BinEdges', -10:0.25:10);
xlabel('zMI Difference (Air-On - Air-Off)');
ylabel('Number of Neurons');
title('C2: Aggregate zMI Differences');

% Plot the distribution of zMI differences for C7
subplot(1, 2, 2);
histogram(zMI_differences_C7, 'BinEdges', -10:0.25:10);
xlabel('zMI Difference (Air-On - Air-Off)');
ylabel('Number of Neurons');
title('C7: Aggregate zMI Differences');

figure(1100);clf;
scatter(zMI_differences_C2,zMI_differences_C7)
%% separate groups based on the difference of zMIs between phases for C2
% Set threshold factor (e.g., 1 for 1 standard deviation)
threshold_factor = 0.5;

% Calculate mean and standard deviation for C2
mean_diff_C2 = mean(zMI_differences_C2, 'omitnan');
std_diff_C2 = std(zMI_differences_C2, 'omitnan');

% Calculate thresholds
threshold_C2_upper = mean_diff_C2 + threshold_factor * std_diff_C2;
threshold_C2_lower = mean_diff_C2 - threshold_factor * std_diff_C2;

% Loop through each animal and calculate the zMI differences for C2 and C7
for ani_idx = 1:numel(zMI_brake_matrices)
    % Get zMI values for the current animal
    zMI_data = zMI_brake_matrices{ani_idx};
    
    % C2 zMI differences (Air-On - Air-Off)
    zMI_diff_C2 = zMI_data(:, 1) - zMI_data(:, 2);
    
    % % C7 zMI differences (Air-On - Air-Off)
    % zMI_diff_C7 = zMI_data(:, 3) - zMI_data(:, 4);
    group_ids_C2_this = NaN(size(zMI_diff_C2));
   group_ids_C2_this(zMI_diff_C2 < threshold_C2_lower) = 1;
   group_ids_C2_this(zMI_diff_C2 > threshold_C2_lower) = 2;
   group_ids_C2_this(isnan(zMI_data(:,1))) = 3;
   group_ids_C2_this(isnan(zMI_data(:,2))) = 4;
   group_ids_C2_this(isnan(zMI_data(:,1)) & isnan(zMI_data(:,2))) = 5;
    group_ids_C2{ani_idx} = group_ids_C2_this;

end

%%
% Initialize matrices to store results for C3, C4, C5, including both Time and Distance zMIs
mean_zMI_C345_group1 = NaN(numel(zMI_no_brake_matrices), 12);  % 12 columns (C3-AOn-T, C3-AOff-T, C4-AOn-T, C4-AOff-T, C5-AOn-T, C5-AOff-T, C3-AOn-D, C3-AOff-D, C4-AOn-D, C4-AOff-D, C5-AOn-D, C5-AOff-D)
mean_zMI_C345_group2 = NaN(numel(zMI_no_brake_matrices), 12);

percent_NaN_group1_C345 = NaN(numel(zMI_no_brake_matrices), 12);
percent_NaN_group2_C345 = NaN(numel(zMI_no_brake_matrices), 12);

% Loop over animals
for ani_idx = 1:numel(zMI_no_brake_matrices)
    % Extract group IDs from C2 (Brake condition) for the current animal
    group_ids_C2_this = group_ids_C2{ani_idx};
    
    % Get zMI values for C3, C4, C5 (No-Brake condition) for the current animal
    zMI_time_C345 = zMI_no_brake_matrices{ani_idx};  % Time encoding data
    zMI_dist_C345 = zMI_no_brake_dist_matrices{ani_idx};  % Distance encoding data
    
    % Group 1 neurons from C2
    group1_idx = (group_ids_C2_this == 1);  % Logical index for group 1 neurons
    group2_idx = (group_ids_C2_this == 2);  % Logical index for group 2 neurons
    
    % Mean zMI for Group 1 in C3, C4, C5 (Time and Distance, Air-On, Air-Off)
    mean_zMI_C345_group1(ani_idx, :) = [mean(zMI_time_C345(group1_idx, :), 1, 'omitnan'), mean(zMI_dist_C345(group1_idx, :), 1, 'omitnan')];
    
    % Mean zMI for Group 2 in C3, C4, C5 (Time and Distance, Air-On, Air-Off)
    mean_zMI_C345_group2(ani_idx, :) = [mean(zMI_time_C345(group2_idx, :), 1, 'omitnan'), mean(zMI_dist_C345(group2_idx, :), 1, 'omitnan')];
    
    % Percent NaN for Group 1 in C3, C4, C5 (Time and Distance, Air-On, Air-Off)
    percent_NaN_group1_C345(ani_idx, :) = [sum(isnan(zMI_time_C345(group1_idx, :)), 1) ./ sum(group1_idx) * 100, sum(isnan(zMI_dist_C345(group1_idx, :)), 1) ./ sum(group1_idx) * 100];
    
    % Percent NaN for Group 2 in C3, C4, C5 (Time and Distance, Air-On, Air-Off)
    percent_NaN_group2_C345(ani_idx, :) = [sum(isnan(zMI_time_C345(group2_idx, :)), 1) ./ sum(group2_idx) * 100, sum(isnan(zMI_dist_C345(group2_idx, :)), 1) ./ sum(group2_idx) * 100];
end

% Now you have:
% - mean_zMI_C345_group1 and mean_zMI_C345_group2 for average zMI across neurons in each group (for C3, C4, C5 with both Time and Distance).
% - percent_NaN_group1_C345 and percent_NaN_group2_C345 for the proportion of NaN values in each group (for C3, C4, C5 with both Time and Distance).

[within,dvn,xlabels] = make_within_table({'G','E','C','P'},[2,2,3,2]);
dataT = make_between_table({[mean_zMI_C345_group1 mean_zMI_C345_group2]},dvn);
raNB_G = RMA(dataT,within,{0.05,{'hsd'}});
%     ra.ranova
print_for_manuscript(raNB_G);

[within,dvn,xlabels] = make_within_table({'G','C','P'},[2,3,2]);
dataT = make_between_table({[mean_zMI_C345_group1(:,7:12) mean_zMI_C345_group2(:,7:12)]},dvn);
raNB_G = RMA(dataT,within,{0.05,{'hsd'}});
%     ra.ranova
print_for_manuscript(raNB_G);