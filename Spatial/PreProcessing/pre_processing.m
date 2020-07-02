% function T = pre_processing
add_to_path
config = get_config;
% clear all
clc
readAgain = 0;
% read excel file if not already read
if ~exist('raw','var') | readAgain
    fileName = 'Data_Info5.xlsx';
    dPath{1} = '\\mohajerani-nas.uleth.ca\storage\homes\samsoon.inayat\Data';
    fileName = fullfile(dPath{1},fileName);
    [num,txt,raw] = xlsread(fileName,10,'A1:N333');
    [f,cName,D] = getFolders;
%     raw = readGoogleSheet;
end
% 
% AD_Thy1_animals = [183224;183227;183228;183329];
AD_Thy1_animals = [183224;183227;183228;183329;1567;1569];
% animals = AD_Thy1_animals;
% Thy1_animals = [173062;173511;173198;174374;173706;183633;183761;183745;183628;183762;1432];
% Thy1_animals = [183628];
animals = AD_Thy1_animals;%
T = getTable(raw,animals,'recording','protocol 15'); 
T = getRecordingFolder(T,D);
T = make_db_and_pdPaths(T,config,2);

T = check_abf_file(T);

T = check_bidiShift(T,0);
% return;
T = check_s2p_run(T);

list = T{:,6};
fileName = 'recording_list.txt';
fid = fopen(fileName,'w');
for ii = 1:length(list)
    try
        thisLine = list{ii};
    catch
        thisLine = list;
    end
    fprintf(fid, '%s\n',thisLine);
end
fclose(fid);

fileName = 'T_16_All.mat';
save(fileName,'T');

