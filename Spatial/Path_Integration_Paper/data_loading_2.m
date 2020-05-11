
% add_to_path
clear all
clc
[f,cName] = getFolders;
T10 = load('T.mat');
selT = T10.T([7 11 13 14 16 18 20 21],:);
ei10 = getData_py(f,selT);
ei10 = loadContextsResponses(ei10,[1 1],[0 1 0]);

selT = T10.T([7 11 13 14 16 18 20 21],:);
ei10 = getData_py(f,selT);
owr = [1,1]; owrp = [1 0 0];
for ii = 1:size(selT,1)
    if ismember(ii,[1:8])
        ei10(ii) = loadContextsResponses(ei10(ii),owr,owrp);
    end
end
disp('Done!');
