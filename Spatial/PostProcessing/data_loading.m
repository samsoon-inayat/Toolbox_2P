% function post_processing
reload = 1;
%% load data table
if reload
    add_to_path
    clear all
    clc
    [f,cName] = getFolders;
    load('T.mat');
    animal_ids = cell2mat(T{:,1});
    ow = 0;
%     eip = getData_py(f,T([7 11 13 15 16 17 18 19 20 21],:));
    indsToLoad = [7 11 13 14 16 18 20 21];
    for ii = 4:length(indsToLoad)
        ii
        eip{ii} = getData_py(f,T(indsToLoad(ii),:));
    end
    disp('Done');
    
    T15 = load('T15.mat');
    animal_ids = cell2mat(T15.T{:,1:2});
    ei15 = getData_py(f,T);
end