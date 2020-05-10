% function post_processing
reload = 1;
%% load data table
if reload
    add_to_path
    clear all
    clc
    [f,cName] = getFolders;
    T_AD = load('T_10_AD.mat');
%     animal_ids = cell2mat(T_AD.T{:,1});
    ow = 0;
    eia = getData_py(f,T_AD.T);
    T_Ctrl = load('T.mat');
    eic = getData_py(f,T_Ctrl.T([17 19 20 21],:));
    disp('Done');
end