% function post_processing
reload = 1;
%% load data table
if reload
%     add_to_path
    clear all
    clc
    [f,cName] = getFolders;
    load('T_15_1_Thy1.mat');
    animal_ids = cell2mat(T{:,1});
    ow = 0;
%     ei = getData(f,T([3 7 11 13 15 16 18 20 21],:));
    ei = getData_py(f,T);
%     eip2 = getData_py(f,T([17 19],:));
    disp('Done');
end

% %% finding all animals
% 
% animal_ids = [cell2mat(T15.T{:,1});cell2mat(T10.T{:,1});cell2mat(T16.T{:,1})];
% uaids = unique(animal_ids);
% length(unique(animal_ids))
% taids = cell2mat(T16.T{:,1});
% inds = [];
% for ii = 1:length(taids)
%     inds(ii) = find(uaids == taids(ii));
% end
% cinds = cellstr(strtrim(num2str(inds')));