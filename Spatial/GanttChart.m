
% %% Define the updated project tasks and their start and end dates
% tasks = {'Planning','Animal Acquisition', 'Surgeries', 'Habituation', 'Training', ...
%          '2P Imaging + EPhys + Behavior', 'Data Analysis', ...
%          'Preparation of Results', 'Preparation of Manuscripts', ...
%          'Presentation at Conferences'};
% 
% % Start and End dates of each task (in datetime format)
% start_dates = datetime({'01-Jul-2025', '01-Aug-2025', '01-Aug-2025', '01-Sep-2025', ...
%                         '01-Oct-2025', '01-Jan-2026', '01-Jan-2027', ...
%                         '01-Mar-2027', '01-Apr-2027', '01-May-2027'}, ...
%                         'InputFormat', 'dd-MMM-yyyy');
% 
% end_dates = datetime({'30-Aug-2025', '30-Jun-2026', '30-Apr-2026', '30-Sep-2025', ...
%                       '31-Dec-2025', '31-Dec-2026', '30-Apr-2027', ...
%                       '31-Mar-2027', '30-Apr-2027', '30-Jun-2027'}, ...
%                       'InputFormat', 'dd-MMM-yyyy');
% 
path = 'G:\My Drive\Research\Grant_Applications\NOT_AG_21_044__PAR_22_094';
path = 'G:\My Drive\Research\Grant_Applications\NOT_AG_24_031_Meta';
% Define the file path and sheet name
file_path = fullfile(path,'TimeLine.xlsx'); % Change to the path of your Excel file
sheet_name = 'Sheet1'; % Specify the sheet name if not the default

% Read the data from the Excel file
data = readtable(file_path, 'Sheet', sheet_name);

% Extract the task names, start dates, and end dates
tasks = data{:, 1}; % Task names in the first column
start_dates = datetime(data{:, 4}, 'InputFormat', 'dd-MMM-yyyy'); % Start dates in the second column
end_dates = datetime(data{:, 7}, 'InputFormat', 'dd-MMM-yyyy'); % End dates in the third column

% % Display the variables
% disp(tasks);
% disp(start_dates);
% disp(end_dates);


% Calculate the duration of each task in days
durations = days(end_dates - start_dates);

% Create a colormap for different task colors
colors = lines(length(tasks)); % Generates distinct colors for each task

% Plot the Gantt chart using rectangles for precise control
% figure(1000);clf

ff = makeFigureRowsCols(107,[3 5 6.9 4],'RowsCols',[1 1],'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.3 0.2],...
    'widthHeightAdjustment',[-350 -400]);

hold on;
for i = 1:length(tasks)
    % Manually plot the bar using a rectangle that adheres to the start and end dates
    rectangle('Position', [datenum(start_dates(i)), i - 0.4, datenum(end_dates(i)) - datenum(start_dates(i)), 0.8], ...
              'FaceColor', colors(i, :), 'EdgeColor', 'none');
end
set(gca, 'YTick', 1:length(tasks), 'YTickLabel', tasks);
xlabel('Absolute Time');
title('Gantt Chart - Project Time Line');
grid on;

% Adjust x-axis to display datetime values
ax = gca;
ax.XLim = [datenum('01-Jun-2025'), datenum('30-Jun-2031')]; % Set limits based on start and end dates
% tick_values = linspace(min(datenum(start_dates)), max(datenum(end_dates)), 3);

datetick('x', 'mmm yy', 'keeplimits'); % Format the x-axis with date labels

% Add task start positions as text
for i = 1:length(tasks)
    text(datenum(start_dates(i)), i, datestr(start_dates(i), 'mmm yy'), ...
         'VerticalAlignment', 'middle', 'HorizontalAlignment', 'right', 'FontSize', 8);
end
x_start = datenum(datetime('01-Dec-2025', 'InputFormat', 'dd-MMM-yyyy')); % X-axis starting point
x_end_line = datenum(datetime('30-Nov-30', 'InputFormat', 'dd-MMM-yyyy')); % Project end line
plot([x_end_line x_end_line],ylim);%, 'End of Project', 'LabelOrientation', 'horizontal', 'LabelVerticalAlignment', 'middle');
plot([x_start x_start],ylim);%, 'End of Project', 'LabelOrientation', 'horizontal', 'LabelVerticalAlignment', 'middle');
% annotation('textbox', [0.4, 0.7, 0.1, 0.1], 'String', 'Project Start', 'FontSize', 8, 'FontWeight', 'bold', 'EdgeColor', 'black', 'BackgroundColor', 'yellow');
hold off;
set(gca,'YDir','reverse','YLim',[0 17]);
xtickangle(15); ytickangle(15);

save_pdf(ff.hf,mData.pdf_folder,sprintf('gantt_chart.pdf'),600);
