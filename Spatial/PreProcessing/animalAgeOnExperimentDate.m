% function animalAgeOnExperimentDate

readAgain = 1;
if readAgain
    clear all
    clc
    fileName = 'Data_Info4.xlsx';
    dPath{1} = '\\mohajerani-nas.uleth.ca\storage\homes\samsoon.inayat\Data';
    fileName = fullfile(dPath{1},fileName);
    [num,txt,raw] = xlsread(fileName,1,'B2:B48');
    animalIDs_list = num;
    [num,txt,raw] = xlsread(fileName,1,'D2:D48');
    DOBs_list = txt;
    T10 = load('T.mat');
    T10_AD = load('T_10_AD.mat');
    T15 = load('T15.mat');
    T16 = load('T16.mat');
    T15_1 = load('T_15_1_Thy1.mat');
    animalIDs = [cell2mat(T10.T{:,1});cell2mat(T15.T{:,1});cell2mat(T16.T{:,1});cell2mat(T10_AD.T{:,1});];
    animalIDs = unique(animalIDs);
    for ii = 1:length(animalIDs)
        ii
        ind = animalIDs_list == animalIDs(ii);
        dobs{ii} = DOBs_list{ind};
    end
end

thisT = T15_1.T;
age = processTable(thisT,animalIDs,dobs); thisT(:,10) = age; thisT.Properties.VariableNames{10} = 'Age';
T15_1.T = thisT;
save('T_15_1_Thy1.mat','-struct','T15_1');

% thisT = T10.T;
% age = processTable(thisT,animalIDs,dobs); thisT(:,10) = age; thisT.Properties.VariableNames{10} = 'Age';
% T10.T = thisT;
% save('T.mat','-struct','T10');
% 
% thisT = T15.T;
% age = processTable(thisT,animalIDs,dobs); thisT(:,10) = age; thisT.Properties.VariableNames{10} = 'Age';
% T15.T = thisT;
% save('T15.mat','-struct','T15');
% 
% thisT = T16.T;
% age = processTable(thisT,animalIDs,dobs); thisT(:,10) = age; thisT.Properties.VariableNames{10} = 'Age';
% T16.T = thisT;
% save('T16.mat','-struct','T16');
% 
% thisT = T10_AD.T;
% age = processTable(thisT,animalIDs,dobs); thisT(:,10) = age; thisT.Properties.VariableNames{10} = 'Age';
% T10_AD.T = thisT;
% save('T_10_AD.mat','-struct','T10_AD');



function age = processTable(T,animalIDs,dobs)
    for ii = 1:size(T,1)
        aID = cell2mat(T{ii,1});
        recordingDate = cell2mat(T{ii,2});
        ind = animalIDs == aID;
        dob = dobs{ind};
        age{ii,1} = 12*years(datetime(recordingDate) - datetime(dob));
    end
end
