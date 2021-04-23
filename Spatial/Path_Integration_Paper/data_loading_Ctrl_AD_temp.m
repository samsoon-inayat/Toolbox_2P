for ii = 1:60
    pause(1);
end
%%
add_to_path
%%
clear all
% clc
%%
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
% \174374\2019-02-12\1_002'
% %%
% ei15_AA(1) = getData_py(f,T15_AD.T(8,:));
% ei15_AA(2) = getData_py(f,T15_AD.T(9,:));
%%
colormaps = load('../MatlabCode/colorblind_colormap.mat');
colormaps.colorblind = flipud(colormaps.colorblind);
mData.colors = mat2cell(colormaps.colorblind,[ones(1,size(colormaps.colorblind,1))]);%{[0 0 0],[0.1 0.7 0.3],'r','b','m','c','g','y'}; % mData.colors = getColors(10,{'w','g'});
mData.axes_font_size = 6; mData.sigColor = [0.54 0.27 0.06];

Uleth_one_drive = 'Z:\homes\brendan.mcallister\2P';
mData.pdf_folder = [Uleth_one_drive '\PDFs']; 
mData.pd_folder = [Uleth_one_drive '\ProcessedDataMatlab'];
disp('Done');
%%
try
    send_email({'samsoon.inayat@uleth.ca'},'Neuroimaging 3 has started loading data for Protocol 10 and 15')
    for ii = 1:size(ET15_C,1)
        ei15_C(ii) = getData_py_old(f,ET15_C(ii,:),0);
    end

    for ii = 1:size(ET15_A,1)
        ei15_A(ii) = getData_py_old(f,ET15_A(ii,:),0);
    end

    ei15_C = loadContextsResponses_ctrl_old(ei15_C,[1 1],[0 0 0]);
    ei15_A = loadContextsResponses_ctrl_old(ei15_A,[1 1],[0 0 0]);
    % ei15_AA = loadContextsResponses(ei15_AA,[1 1],[0 0 0]);
    % ei15_A(3) = ei15_AA(2);


    % for loading behavior and 2p data
    for ii = 1:2%size(ET10_C,1)
        ei10_C(ii) = getData_py_old(f,ET10_C(ii,:),0);
    end

    for ii = 1:2%size(ET10_A,1)
        ei10_A(ii) = getData_py_old(f,ET10_A(ii,:),0);
    end

    ei10_C = loadContextsResponses_ctrl_old(ei10_C,[1 1],[0 0 0]);
    ei10_A = loadContextsResponses_ctrl_old(ei10_A,[1 1],[0 0 0]);
    training_data_C = behaviorProcessor;
    training_data_A = behaviorProcessor_AD;

    parameter_matrices_ctrl('calculate','10_C',ei10_C);
    parameter_matrices_ctrl('calculate','15_C',ei15_C);
    parameter_matrices_ctrl('calculate','10_A',ei10_A);
    parameter_matrices_ctrl('calculate','15_A',ei15_A);
%     
%     parameter_matrices('calculate','10_C',ei10_C);
%     parameter_matrices('calculate','15_C',ei15_C);
%     parameter_matrices('calculate','10_A',ei10_A);
%     parameter_matrices('calculate','15_A',ei15_A);
    
    send_email({'samsoon.inayat@uleth.ca','brendan.mcallister@uleth.ca'},'Complete - Loading data for Protocol 10 and 15')
catch
    send_email({'samsoon.inayat@uleth.ca'},'Error occurred while loading data')
end
%%
for ii = 1:size(ET15,1)
    ei15(ii) = getData_py(f,ET15(ii,:));
end


%%
ei10 = loadContextsResponses_1(ei10,[1 1],[1 -1 1]);
ei15 = loadContextsResponses(ei15,[1 1],[0 0 0]);
ei16 = loadContextsResponses(ei16,[1 1],[0 0 0]);

% parameter_matrices('calculate','10_C',ei10_C);
% parameter_matrices('calculate','15_C',ei15_C);
% parameter_matrices('calculate','16',ei16);
disp('All Done!');
%%


ei10 = loadContextsResponses(ei10,[0 0],[-1 -1 -1]);
ei10 = loadContextsResponses(ei10,[1 1],[0 0 0]);
% ei10 = loadContextsResponses(ei10,1,[1 1 1]);

%% for Sam-WS
owr = [1,1]; owrp = [0 0 0];
for ii = 1:length(selRecs)
    if ismember(ii,[8])
        ei10(ii) = loadContextsResponses(ei10(ii),owr,owrp);
    end
end
disp('Done!');
%%
glms = do_glm(ei10,0);
glmsI = do_glmI(ei10,0);
disp('Done!');

%% for neuroimaging computer
owr = [1,1]; owrp = [0 1 1];
for ii = 1:size(selT,1)
    if ismember(ii,[1:8])
        ei10(ii) = loadContextsResponses(ei10(ii),owr,owrp);
    end
end
disp('Done!');

owr = [1,1]; owrp = [1 0 0];
for ii = 1:size(selT,1)
    if ismember(ii,[1:8])
        ei10(ii) = loadContextsResponses(ei10(ii),owr,owrp);
    end
end
disp('Done!');
%%
glms = do_glm(ei10(9),1);
glmsI = do_glmI(ei10(9),1);
disp('Done!');

%%

T15 = load('T15.mat');
ei15 = getData_py(f,T15.T([8 2 6 4 10],:));


%%
% this is just to load behavior

for ii = 1:size(T10.T,1)
    if ismember(ii,[1:size(T10.T,1)]) % select which data to load in the second argument
        eiB(ii) = getBehavior(f,T10.T(ii,:));
    end
end
disp('Done!');

%%
% This is just to visualize behavior graphs
inds = [];
for ii = 1:size(T10.T,1)
    if ismember(ii,[1:size(T10.T,1)]) % select which data to load in the second argument
        if isempty(eiB{ii})
            continue;
        end
         behaviorPlot(eiB(ii))
         ii
         key = getkey;
         if key == 27 % esc
             break;
         end
         if key == 105 %i
             inds = [inds ii];
         end
    end
end
inds
disp('Done!');

