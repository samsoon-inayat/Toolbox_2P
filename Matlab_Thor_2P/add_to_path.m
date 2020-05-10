restoredefaultpath

[~, name] = system('hostname');

if strcmp(strtrim(name),'Sam-WS')
%     scriptPath = 'G:\OneDrives\OneDrive\OneDrive_U_Leth\Google Drive';
    scriptPath = 'G:\OneDrives\OneDrive - University of Lethbridge';
end

if strcmp(strtrim(name),'imaging3-psy16')
    scriptPath = 'E:\Users\samsoon.inayat\OneDrive - University of Lethbridge';
end

if strcmp(strtrim(name),'imaging4-psy16')
    scriptPath = 'E:\Users\samsoon.inayat\OneDrive - University of Lethbridge';
end

if strcmp(strtrim(name),'imaging1-psy16')
    scriptPath = 'E:\Users\samsoon.inayat\OneDrive - University of Lethbridge';
end

% addpath(scriptPath);
matlabCodeFolder = fullfile(scriptPath,'MatlabCode');addpath(matlabCodeFolder);
temp = fullfile(matlabCodeFolder,'Common');addpath(temp);
temp = fullfile(matlabCodeFolder,'Common\append_pdfs');addpath(temp);
temp = fullfile(matlabCodeFolder,'Thor2P');addpath(temp);
temp = fullfile(matlabCodeFolder,'FigureFunctions');addpath(temp);
temp = fullfile(matlabCodeFolder,'ca_source_extraction-master');p = genpath(temp);addpath(p);

temp = fullfile(matlabCodeFolder,'Suite2P-master');p = genpath(temp);addpath(p);

temp = pwd; p = genpath(temp);addpath(p);
remove_from_path
% addToPath(matlabCodeFolder);
% rmpath(scriptPath);


