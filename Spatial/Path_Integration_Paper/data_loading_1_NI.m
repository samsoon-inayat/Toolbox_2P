
% add_to_path
clear all
clc
mData.colors = {[0 0 0],[0.1 0.7 0.3],'r','b','m','c','g','y'};
% mData.colors = getColors(10,{'w','g'});
mData.axes_font_size = 6;
mData.sigColor = [0.54 0.27 0.06];
mData.pdf_folder = fullfile(pwd,'PDFs');
disp('data extracted');

[f,cName] = getFolders;
T10 = load('T_10_All.mat');
T15 = load('T_15_All.mat');
T16 = load('T_16_All.mat');

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


%%
% for loading behavior and 2p data
selRecs = [4     8    12    15    16    17    18    19    20    21    22    24    25];
% selRecs = selRecs([2 3 4 5 7 9 11 12 13]);
T10.T{selRecs,1};
ind = 1;
for ii = 1:size(T10.T,1)
    if ismember(ii,selRecs) % select which data to load in the second argument
        ei10(ind) = getData_py(f,T10.T(ii,:));
        ind = ind + 1;
    end
end
disp('Done');

%% for loading behavior and 2p data
selRecs = [1:9 12 13 16];
% selRecs = [10 11 14 15];
T15.T{selRecs,4};
ind = 1;
for ii = 1:size(T15.T,1)
    disp(sprintf('%d of %d',ii,size(T15.T,1)));
    if ismember(ii,selRecs) % select which data to load in the second argument
        ei15(ind) = getData_py(f,T15.T(ii,:));
        ind = ind + 1;
    end
end
disp('Done');

%% for loading behavior and 2p data
selRecs = [1 2 4 5];
% selRecs = [10 11 14 15];
T16.T{selRecs,4};
ind = 1;
for ii = 1:size(T16.T,1)
    disp(sprintf('%d of %d',ii,size(T16.T,1)));
    if ismember(ii,selRecs) % select which data to load in the second argument
        ei16(ind) = getData_py(f,T16.T(ii,:));
        ind = ind + 1;
    end
end
disp('Done');

%%
ei10(8) = getData_py(f,T10.T(selRecs(8),:));

%%
% ei10 = getData_py(f,selT);
ei10 = loadContextsResponses(ei10,[1 1],[0 0 0]);

disp('Done');

%%
% ei10 = getData_py(f,selT);
ei16 = loadContextsResponses(ei16,[1 1],[0 0 0]);

disp('Done');

%%
parameter_matrices('calculate',ei10);
disp('Done!');
%%
training_data = behaviorProcessor;

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

