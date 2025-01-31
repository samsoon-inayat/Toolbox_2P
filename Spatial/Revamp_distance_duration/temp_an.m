
an = 1;
tudata = udataT{an}
field_names = fieldnames(tudata);
for ii = 1:length(field_names)
    varname = field_names{ii};
    cmdTxt = sprintf('%s = tudata.%s;',varname,varname);
    eval(cmdTxt);
end
figure(1000);clf;
plot(ts,animal_motion);hold on
plot(ts,frf);
% plot(ts,bnb);
% plot(ts,C1);
% plot(ts,C2);
% plot(ts,C3);
% plot(ts,C4);
% plot(ts,C5);
% plot(ts,C6);
% plot(ts,C7);

%%
% Define activation threshold (mean + 2*std across all time for each neuron)
threshold = mean(firing_rate, 2) + 2 * std(firing_rate, 0, 2);

% Logical indices for brake and no-brake conditions
brake_idx = bnb == 1;   % Brake conditions
nobrake_idx = bnb == 0; % No-brake conditions

% Logical indices for air-on and air-off phases
airon_idx = air == 1;
airoff_idx = air == 0;

% Preallocate results
results = struct('brake_airon', [], 'brake_airoff', [], 'nobrake_airon', [], 'nobrake_airoff', []);

% Calculate percentage of activated cells
% Brake + Air-On
active_brake_airon = mean(firing_rate(:, brake_idx & airon_idx) > threshold, 2);
results.brake_airon = 100 * sum(active_brake_airon > 0) / size(firing_rate, 1);

% Brake + Air-Off
active_brake_airoff = mean(firing_rate(:, brake_idx & airoff_idx) > threshold, 2);
results.brake_airoff = 100 * sum(active_brake_airoff > 0) / size(firing_rate, 1);

% No-Brake + Air-On
active_nobrake_airon = mean(firing_rate(:, nobrake_idx & airon_idx) > threshold, 2);
results.nobrake_airon = 100 * sum(active_nobrake_airon > 0) / size(firing_rate, 1);

% No-Brake + Air-Off
active_nobrake_airoff = mean(firing_rate(:, nobrake_idx & airoff_idx) > threshold, 2);
results.nobrake_airoff = 100 * sum(active_nobrake_airoff > 0) / size(firing_rate, 1);

% Display results
disp('Percentage of Activated Cells:');
disp(results);

% Optional: Bar Plot of Results
figure;
bar([results.brake_airon, results.brake_airoff, results.nobrake_airon, results.nobrake_airoff]);
xticks(1:4);
xticklabels({'Brake + Air-On', 'Brake + Air-Off', 'No-Brake + Air-On', 'No-Brake + Air-Off'});
ylabel('Percentage of Activated Cells (%)');
title('Percentage of Activated Cells in Different Conditions');
