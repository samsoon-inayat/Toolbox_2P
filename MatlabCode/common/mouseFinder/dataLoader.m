% Written by Sam Inayat on August 07, 2018
% Opens the excel file containing list of mice and by choosing different
% parameters such as age, genotype, gender, one can find desired list of
% mice

clear all
clc
% function miceFinder
path = '\\mohajerani-nas.uleth.ca\storage\mice_booking';
% fileName = 'Thy1-GCaMP.xlsx';
% fileName = 'Thy1-GCaMP.xlsx';
fileName = 'APP-Thy1GCaMP.xlsx';

fileName = fullfile(path,fileName);

lastRow = 264;
lastColumn = 'P';

xlsreadRange = sprintf('A2:%s%d',lastColumn,lastRow);
[num,txt,raw] = xlsread(fileName,1,xlsreadRange);
