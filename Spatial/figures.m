function figures
SamSamSamSamSamSamSam % put here to create an error o that this file never runs ... only for sequence definition
% I will never name a variable like that

%%
dataLoader1
%% initial data processing (once run ... only run again if need processing again)
saving_mean_raster_fits

%% figure 1B
behaviorProcessor_1

%% figure 1C
plotSpeedTrialsImage
plotTimeToOnsetOfMovementAfterAirPuff
plotAverageSpeedConditions

%% figure 1D
speedFigure



%% figure 1F
tracePlots_1

%% figure 1G
tracePlots_1


%% finding place cells
figure_place_cells_vs_other_cells
pdf = {'figure_place_cells.pdf'};

%% place field properties (distributions, MI, Rsq, PW, PC and more)
figure_place_field_properties
pdf = {'figure_pf_widths_dist.pdf','figure_pf_centers_dist.pdf','figure_scatter_centers_widths.pdf','figure_widths_vs_centers.pdf'};
pdf = {'figure_MI_dist.pdf','figure_RS_dist.pdf'};

%% place field dynamics
figure_place_field_dynamics_trials

figure_place_field_dynamics_contexts
type1 = 'cells of context 1 disrupted in context 2';
type2 = 'cells of context 1 that did not get disrupted in context 2';