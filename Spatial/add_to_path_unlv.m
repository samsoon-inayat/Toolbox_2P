function add_to_path
currentDir = pwd;
strpos = findstr(currentDir,'GitHub');
github_dir = currentDir(1:(strpos+length('GitHu')));
interDir = fullfile(github_dir,'Common');
mainDir = pwd;
% if ~exist(fullfile(mainDir,'Downloaded'),'dir')
    % mainDir = '\\mohajerani-nas.uleth.ca\storage\homes\brendan.mcallister\2P\Analysis';
    mainDir = 'E:\GoogleDrive\InayatSamsoon\ULETH\Analysis';
% end
cd(currentDir);
% addpath(scriptPath);

temp = fullfile(interDir,'Matlab');addpath(temp);
temp = fullfile(interDir,'Matlab','Common');addpath(temp);
temp = fullfile(interDir,'Matlab','Downloaded');addpath(temp);
temp = fullfile(mainDir,'Downloaded','append_pdfs');addpath(temp);
temp = fullfile(interDir,'Matlab','FigureFunctions');addpath(temp);
temp = fullfile(mainDir,'Downloaded','spatial-information-metrics');addpath(temp);
temp = fullfile(interDir,'Matlab','Fitting');addpath(temp);
temp = fullfile(interDir,'Matlab','Statistical_Package');addpath(temp);
temp = fullfile(mainDir,'Downloaded','npy-matlab');addpath(temp);
temp = fullfile(mainDir,'Downloaded');addpath(temp);
temp = fullfile(mainDir,'Downloaded','CNMF_E');p = genpath(temp);addpath(p);
temp = fullfile(mainDir,'Downloaded','Suite2P');p = genpath(temp);addpath(p);

% temp = pwd; p = genpath(temp);addpath(p);
% remove_from_path

folders_to_add = {'Misc','ProcessData','ProcessedDataIO','Thor2P','PreProcessing','contextsDefineAndProcess','Visualization'};

for ii = 1:length(folders_to_add)
    p = genpath(fullfile(pwd,folders_to_add{ii}));
    addpath(p);
end


