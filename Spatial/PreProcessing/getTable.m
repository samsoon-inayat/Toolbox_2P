function T = getTable(raw,animals,TR,P)

T = [];
for ii = 1:length(animals) % for all the animals in the animal list make a table for the selected protocol or recording /training
    sel_id = animals(ii);
    TT = getTable_1(raw,sel_id,TR,P); 
    T = [T;TT];
end


function T = getTable_1(raw,sel_id,TR,P)

% columns in excel file
IDs = getSelectedCol(raw,'A');
RecordingDate = getSelectedCol(raw,'B');
TrainingOrRecording = getSelectedCol(raw,'E');
Protocol = getSelectedCol(raw,'I');
Quality = getSelectedCol(raw,'J');
PerformanceNotes = getSelectedCol(raw,'N');

selection = [];



% specify min and max age
ii = 1;
ids = getIDs(raw,'A');
ids = ids';
selection{ii} = find(ids == sel_id); ii = ii + 1;
if ~isempty(TR)
    IndexC = strfind(lower(TrainingOrRecording), TR);
    selection{ii} = find(not(cellfun('isempty', IndexC))); ii = ii + 1;
end
if ~isempty(P)
    IndexC = strfind(lower(Protocol), P);
    selection{ii} = find(not(cellfun('isempty', IndexC))); ii = ii + 1;
end

% common = intersect(selection{1},selection{2});
common = selection{1};
if length(selection) > 1
    for ii = 2:length(selection)
        common = intersect(common,selection{ii});
    end
end
bad = [];
for ii = 1:length(common)
    if ~isempty(strfind(lower(PerformanceNotes{common(ii)}),'euthanize')) || ~isempty(strfind(lower(PerformanceNotes{common(ii)}),'bad'))...
            || ~isempty(strfind(lower(Quality{common(ii)}),'no data')) || ~isempty(strfind(lower(Quality{common(ii)}),'do not'))
        bad = [bad ii];
    end
end
common(bad) = [];
% 

varNames = {'ID','RecordingDate','TrainingOrRecording','Protocol','PerformanceNotes'};
T = table(ids(common),RecordingDate(common),TrainingOrRecording(common),Protocol(common),PerformanceNotes(common),'VariableNames',varNames);

