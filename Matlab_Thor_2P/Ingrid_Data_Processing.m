% Ingrid's Data Processing for checking my spatial information calculator
% code

clear all
clc
rootFolder = '\\mohajerani-nas.uleth.ca\storage2\homes\samsoon.inayat\Data\RSC037 11_10';

files = {'behavior.mat','deconv.mat','PC_list.mat','timecourses'};
for ii = 1:length(files)
    fileName = files{ii};%makeName(files{ii},rootFolder);
    load(fileName)
end

% onsets = 

for ii = 1:size(deconv,2)
    
    
end