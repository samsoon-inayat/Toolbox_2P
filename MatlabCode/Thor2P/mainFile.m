% flow of the code

% 1. Take the raw data and find a good set of frames where there is no
% motion ... then find the average of those frames and use the average to
% register other frames
% 2. After motion correcting frames, down sample them so that they can be
% used for cross correlation calculations
% 3. Also find the average image from motion corrected frames

clear all
clc

%  listAllData('F:\Sam\Data');

% stepsToExecuteNames = {'MotionCorrection';'AverageImage';'DownSampling';'IdentifyROIs';'ROIsSignals'};

% 1. Select the data that you want to process
% dataFolder = 'F:\Sam\Data\160513\CoverSlip940\06_10_2016\EPhys';
% % dataFolder = 'F:\Sam\Data\Animal_160513\CoverSlip940\06_10_2016\EPhys';
% % dataFolder = 'G:\Data\160511\06_13_2016\EPhys';
% % dataFolder = 'F:\Sam\Data\Animal_153465\05_18_2016\freeRun003';
% dataFolder = 'F:\Sam\Data\160558\CoverSlip\06_09_2016\EPhys2P001';
% dataFolder = 'F:\Sam\Data\160511\06_13_2016\EPhys';
dataFolder = 'F:\Sam\Data\160558\07_08_2016\EPhys2P001';
ei = thorGetExperimentInfo(dataFolder);
% 
% % movieMaker(ei,[1:500],'tif');
% 
raw2tif(ei,[1 1000]);
makeRegisterationImage(ei);
motionCorrection(ei);
findAverageImage(ei);
% % downSampleRawFile(ei);
saveROIsSignal(ei);
% 
processABF(ei);
plotBehavior(ei);
plotROIs_all_in_one_figure(ei)
plotROIs(ei)