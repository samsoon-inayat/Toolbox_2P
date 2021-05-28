function temporary_contexts_translator

[f,cName] = getFolders;
T10 = load('T_10_All.mat');
T15 = load('T_15_All.mat');
selRecs10 = [4     8    12    15    16    17    18    19    20    21    22    24    25];
% selRecs15 = [1:9 12 13 16];
selRecs15 = [1 2 4 6 8 12 16];
ET10_C = T10.T(selRecs10,:); ET15_C = T15.T(selRecs15,:); 
ET10_C = ET10_C([6 7 9 11 12],:);

T10_AD = load('T_10_All_AD.mat');
T15_AD = load('T_15_All_AD.mat');
ET10_A = T10_AD.T(2:6,:); ET15_A = T15_AD.T([1 4 9 10 13 14],:); 
% clear('selRecs10','selRecs15','T10','T15','T10_AD','T15_AD','cName')
disp('done')

newTC = evalin('base','T_C');
newTA = evalin('base','T_A');


% tei{ii}.plane{pp}.folder = fullfile(plane{pp},'post_suite2p_matlab');

% for ii = 1:size(newTC,1)
%     oldcfile = fullfile(cell2mat(ET10_C{ii,7}),'suite2p\plane0','post_suite2p_matlab','contexts.mat');
%     temp = load(oldcfile);
%     newfile = fullfile(cell2mat(ET10_C{ii,6}),'contexts.mat');
%     if ~copyfile(oldcfile,newfile);
%         n = 0;
%     end
% end

for ii = 1:size(newTA,1)
    oldcfile = fullfile(cell2mat(ET10_A{ii,7}),'suite2p\plane0','post_suite2p_matlab','contexts.mat');
    temp = load(oldcfile);
    newfile = fullfile(cell2mat(ET10_A{ii,6}),'contexts.mat');
    if ~copyfile(oldcfile,newfile);
        n = 0;
    end
end



