% function post_processing
reload = 1;
%% load data table
if reload
%     add_to_path
    clear all
    clc
    [f,cName] = getFolders;
    load('T16.mat');
    animal_ids = cell2mat(T{:,1});
    ow = 0;
%     ei = getData(f,T([3 7 11 13 15 16 18 20 21],:));
    ei = getData_py(f,T);
%     eip2 = getData_py(f,T([17 19],:));
    disp('Done');
end