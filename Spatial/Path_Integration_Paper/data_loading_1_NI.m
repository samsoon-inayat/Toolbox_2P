
% add_to_path
clear all
clc
[f,cName] = getFolders;
T10 = load('T.mat');
selT = T10.T([7 11 13 14 16 18 20 21 22],:);
for ii = 1:size(selT,1)
    if ismember(ii,[1:9]) % select which data to load in the second argument
        ei10(ii) = getData_py(f,selT(ii,:));
    end
end
% ei10 = getData_py(f,selT);
ei10 = loadContextsResponses(ei10,[1 1],[0 0 0]);

disp('Done');
%%
training_data = behaviorProcessor;

ei10 = loadContextsResponses(ei10,[0 0],[-1 -1 -1]);
ei10 = loadContextsResponses(ei10,[1 1],[0 0 0]);
% ei10 = loadContextsResponses(ei10,1,[1 1 1]);

%% for Sam-WS
owr = [1,1]; owrp = [0 0 0];
for ii = 1:size(selT,1)
    if ismember(ii,[9])
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
mData.colors = {[0 0 0],[0.1 0.7 0.3],'r','b','m','c','g','y'};
% mData.colors = getColors(10,{'w','g'});
mData.axes_font_size = 6;
mData.sigColor = [0.54 0.27 0.06];
mData.pdf_folder = fullfile(pwd,'PDFs');
disp('data extracted');