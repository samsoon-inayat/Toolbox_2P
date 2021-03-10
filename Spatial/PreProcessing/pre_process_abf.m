% [f,cName] = getFolders;
% load('T_15_1_Thy1.mat');

process_abf(T15_c(sel15,:));


ei15_AD = getData_py(f,T15_AD.T(13,:));
processContextDefinitions(ei15_AD);