for ii = 1:1000
    pause(1);
end

add_to_path
%%
% add_to_path
% clear all
% clc
[f,cName] = getFolders;
T10 = load('T_10_All.mat');
T15 = load('T_15_All.mat');
T16 = load('T_16_All.mat');
selRecs10 = [4     8    12    15    16    17    18    19    20    21    22    24    25];
selRecs15 = [1:9 12 13 16];
selRecs16 = [1 2 4 5];
ET10 = T10.T(selRecs10,:); ET15 = T15.T(selRecs15,:); ET16 = T16.T(selRecs16,:);
disp('Done');

%%
colormaps = load('../MatlabCode/colorblind_colormap.mat');
colormaps.colorblind = flipud(colormaps.colorblind);
mData.colors = mat2cell(colormaps.colorblind,[ones(1,size(colormaps.colorblind,1))]);%{[0 0 0],[0.1 0.7 0.3],'r','b','m','c','g','y'}; % mData.colors = getColors(10,{'w','g'});
mData.axes_font_size = 6; mData.sigColor = [0.54 0.27 0.06]; mData.pdf_folder = fullfile(pwd,'PDFs'); 
mData.selAnimals10_d1 = [2 3 4 5 11];%[1:5 7 9 11:13]; 
mData.selAnimals10_af15 = [7 9 11 12];%[1:5 7 9 11:13]; 
mData.selAnimals10 = mData.selAnimals10_d1;
mData.selAnimals10 = mData.selAnimals10_af15;
mData.selAnimals10 = [mData.selAnimals10_af15 mData.selAnimals10_d1];
mData.selAnimals15_d1 = [8 2 6 4 10];%[1 2 4 6 8 10 12]; 
mData.selAnimals15_d2 = mData.selAnimals15_d1 + 1;%[1 2 4 6 8 10 12]; 
mData.selAnimals15 = mData.selAnimals15_d1;
mData.selAnimals15 = mData.selAnimals15_d2;
mData.selAnimals15 = [mData.selAnimals15_d1 mData.selAnimals15_d2];
mData.selAnimals16 = [1 2 3 4];
disp('Done');

%%
% for loading behavior and 2p data
for ii = 1:size(ET10,1)
    ei10(ii) = getData_py(f,ET10(ii,:));
end

for ii = 1:size(ET15,1)
    ei15(ii) = getData_py(f,ET15(ii,:));
end

for ii = 1:size(ET15_AD,1)
    ei15(ii) = getData_py(f,ET15(ii,:));
end


for ii = 1:size(ET16,1)
    ei16(ii) = getData_py(f,ET16(ii,:));
end

ei10 = loadContextsResponses(ei10,[1 1],[0 0 0]);
ei15 = loadContextsResponses(ei15,[1 1],[0 0 0]);
ei16 = loadContextsResponses(ei16,[1 1],[0 0 0]);

parameter_matrices('calculate','10',ei10);
parameter_matrices('calculate','15',ei15);
parameter_matrices('calculate','16',ei16);
disp('All Done!');
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

