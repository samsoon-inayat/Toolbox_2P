function add_to_path
currentDir = pwd;
cd('..');
interDir = pwd;
cd('..');
mainDir = pwd;
cd(currentDir);
% addpath(scriptPath);
temp = fullfile(interDir,'MatlabCode');addpath(temp);
temp = fullfile(interDir,'MatlabCode','Common');addpath(temp);
temp = fullfile(mainDir,'Downloaded','append_pdfs');addpath(temp);
temp = fullfile(interDir,'MatlabCode','Thor2P');addpath(temp);
temp = fullfile(interDir,'MatlabCode','FigureFunctions');addpath(temp);
temp = fullfile(mainDir,'Downloaded','spatial-information-metrics-master');addpath(temp);
temp = fullfile(mainDir,'Downloaded','npy-matlab-master');addpath(temp);
temp = fullfile(mainDir,'Downloaded');addpath(temp);
temp = fullfile(mainDir,'Downloaded','ca_source_extraction-master');p = genpath(temp);addpath(p);


temp = fullfile(mainDir,'Downloaded','Suite2P-master');p = genpath(temp);addpath(p);

temp = pwd; p = genpath(temp);addpath(p);
remove_from_path



