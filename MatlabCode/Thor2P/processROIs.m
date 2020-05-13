clear all
clc

% dataFolder = 'F:\Sam\Data\Animal_153465\05_17_2016\freeRun';
dataFolder = 'F:\Sam\Data\Animal_160513\CoverSlip940\05_29_2016\freeRun';
expInfo = thorGetExperimentInfo(dataFolder);
ROIs = getROIPixels(expInfo);

display('done');