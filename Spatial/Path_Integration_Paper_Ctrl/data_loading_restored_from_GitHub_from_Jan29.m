for ii = 1:60
    pause(1);
end
clear all
%%
% add_to_path

% clc
[f,cName] = getFolders;
T10 = load('T_10_All.mat');
T15 = load('T_15_All.mat');
selRecs10 = [4   8    12    15    16    17    18    19    20    21    22    24    25];
% selRecs15 = [1:9 12 13 16];
% selRecs15 = [1 2 4 6 8 12 16];
ET10_C = T10.T(selRecs10,:); %ET15_C = T15.T(selRecs15,:); 
ET10_CC = ET10_C([6 7 9 11 12],:);
ET10_CC1 = ET10_C([4 5 13],:);
ET10_CC1{1,7} = {'\\mohajerani-nas.uleth.ca\storage\homes\brendan.mcallister\2P\Processed_Data\OLD\174374\2019-02-12\1_002'};
ET10_CC1{2,7} = {'\\mohajerani-nas.uleth.ca\storage\homes\brendan.mcallister\2P\Processed_Data\OLD\173706\2019-02-08\1'};
ET10_CC1{3,7} = {'\\mohajerani-nas.uleth.ca\storage\homes\brendan.mcallister\2P\Processed_Data\OLD\001432\2020-05-02\1_001'};
ET10_CD = ET10_C([1 2 3],:);
ET10_Comb = [ET10_CD;ET10_CC];

% T10_AD = load('T_10_All_AD.mat');
% T15_AD = load('T_15_All_AD.mat');
% ET10_A = T10_AD.T(2:6,:); ET15_A = T15_AD.T([1 4 9 10 13 14],:); 
% clear('selRecs10','selRecs15','T10','T15','T10_AD','T15_AD','cName')
disp('done')
% \174374\2019-02-12\1_002'
% %%
% ei15_AA(1) = getData_py(f,T15_AD.T(8,:));
% ei15_AA(2) = getData_py(f,T15_AD.T(9,:));
%%
colormaps = load('../MatlabCode/colorblind_colormap.mat');
colormaps.colorblind = flipud(colormaps.colorblind);
% mData.colors = mat2cell(colormaps.colorblind,[ones(1,size(colormaps.colorblind,1))]);%
mData.colors = {[0 0 0],[0.1 0.7 0.3],'r','b','m','c','g','y'}; % mData.colors = getColors(10,{'w','g'});
mData.axes_font_size = 6; mData.sigColor = [0.54 0.27 0.06]; 
Uleth_one_drive = 'G:\OneDrives\OneDrive - University of Lethbridge';
Uleth_one_drive = 'E:\Users\samsoon.inayat\OneDrive - University of Lethbridge';
mData.pdf_folder = [Uleth_one_drive '\PDFs']; 
mData.pd_folder = [Uleth_one_drive '\ProcessedData\Matlab'];
disp('Done');
%%
try
    send_email({'samsoon.inayat@uleth.ca'},'Neuroimaging 3 has started loading data for Protocol 10 and 15')
%     for ii = 1:size(ET15_C,1)
%         ei15_C(ii) = getData_py(f,ET15_C(ii,:));
%     end
% 
% %     for ii = 1:size(ET15_A,1)
% %         ei15_A(ii) = getData_py(f,ET15_A(ii,:));
% %     end
% 
%     ei15_C = loadContextsResponses_ctrl(ei15_C,[1 1],[0 0 0]);
%     ei15_A = loadContextsResponses_ctrl(ei15_A,[1 1],[0 0 0]);
    % ei15_AA = loadContextsResponses_ctrl(ei15_AA,[1 1],[0 0 0]);
    % ei15_A(3) = ei15_AA(2);


    % for loading behavior and 2p data
    for ii = 1:size(ET10_CD,1)
        ei10_C(ii) = getData_py(f,ET10_CD(ii,:),0);
    end
    
%     for ii = 1:size(ET10_CC1,1)
%         ei10_C1(ii) = getData_py(f,ET10_CC1(ii,:));
%     end

    for ii = 1:size(ET10_CC,1)
        ei10_A(ii) = getData_py(f,ET10_CC(ii,:),0);
    end

    ei10_C = loadContextsResponses_ctrl_old(ei10_C,[1 1],[0 0 0]);
%     ei10_C1 = loadContextsResponses_ctrl(ei10_C1,[1 1],[0 0 0]);
    ei10_A = loadContextsResponses_ctrl_old(ei10_A,[1 1],[0 0 0]);
    ei_comb = [ei10_C ei10_A];
    ei_150 = ei10_A(2:5);
    training_data_C1 = behaviorProcessor;
%     1432      173062      173198      173511      173706      174374      183628      183633      183745      183761      183762
    training_data_C1.belt_lengths = [150 142 142 142 142 142 150 150 150 150 150]';
    weight_day1 = [52.1 46.7 47.3 39.9 33.7 30.9 26.8 34.5 34.6 31.6 35.8]';
    weight_day2 = [52.8 46.4 47.8 40.2 33.3 30.4 26.6 34.2 34.5 31.3 35.4]';
    weight_day3 = [51.6 45.9 47.9 39.5 33.3 30.5 26.7 33.7 33.6 30.7 35.3]';
    weight_day4 = [50.4 45.9 46.4 37.8 31.8 29.7 27.1 33.4 36.1 31.5 35.8]';
    training_data_C1.weight = [weight_day1 weight_day2 weight_day3 weight_day4];
    training_data_C1.DOB = {'2018-09-27','2017-08-10','2017-08-28','2017-09-18','2017-10-15','2017-12-09',...
        '2018-09-27','2018-09-27','2018-10-03','2018-10-11','2018-10-11'};
%     training_data_A = behaviorProcessor_AD;
    for rr = 1:size(training_data_C1.training_dates,1)
        dob = datetime(training_data_C1.DOB{rr},'InputFormat','yyyy-MM-dd');
        for cc = 1:size(training_data_C1.training_dates,2)
            this_date = training_data_C1.training_dates(rr,cc); this_date = this_date{1};
            if isempty(this_date)
                ageAn(rr,cc) = -1;
            else
                dv = datevec(this_date-dob);
                ageAn(rr,cc) = 12*dv(1) + dv(2) + dv(3)/30;
            end
        end
    end
    animal_id_C = [173062,173511,173198];
    date_of_rec_C = {'2018-07-11','2018-07-11','2018-07-11'};
    date_of_surg_C = {'2018-04-13','2018-04-30','2018-04-27'};
    animal_id_A = [173706,183761,183745,183628,183762];
    date_of_rec_A = {'2019-02-14','2019-06-19','2019-06-19','2019-06-21','2019-06-24'};
    date_of_surg_A = {'2018-10-15','2019-03-29','2019-04-30','2019-04-03','2019-04-26'};
    for rr = 1:length(animal_id_C)
        ind = find(training_data_C1.animalIDs == animal_id_C(rr));
        dob = datetime(training_data_C1.DOB{ind},'InputFormat','yyyy-MM-dd');
        this_date = datetime(date_of_rec_C{rr},'InputFormat','yyyy-MM-dd');
        dv = datevec(this_date-dob);
        ageAn_rec_C(rr) = 12*dv(1) + dv(2) + dv(3)/30;
        start_date = datetime(training_data_C1.training_dates{ind,1},'InputFormat','yyyy-MM-dd');
        dv = datevec(this_date-start_date);
        rec_day_C(rr) = 12*dv(1) + dv(2) + dv(3)/30;
        surg_date = datetime(date_of_surg_C{rr},'InputFormat','yyyy-MM-dd');
        dv = datevec(start_date-surg_date);
        postsurg_day_C(rr) = 12*dv(1) + dv(2) + dv(3)/30;
    end
    for rr = 1:length(animal_id_A)
        ind = find(training_data_C1.animalIDs == animal_id_A(rr));
        dob = datetime(training_data_C1.DOB{ind},'InputFormat','yyyy-MM-dd');
        this_date = datetime(date_of_rec_A{rr},'InputFormat','yyyy-MM-dd');
        dv = datevec(this_date-dob);
        ageAn_rec_A(rr) = 12*dv(1) + dv(2) + dv(3)/30;
        start_date = datetime(training_data_C1.training_dates{ind,1},'InputFormat','yyyy-MM-dd');
        dv = datevec(this_date-start_date);
        rec_day_A(rr) = 12*dv(1) + dv(2) + dv(3)/30;
        surg_date = datetime(date_of_surg_A{rr},'InputFormat','yyyy-MM-dd');
        dv = datevec(start_date-surg_date);
        postsurg_day_A(rr) = 12*dv(1) + dv(2) + dv(3)/30;
    end

    parameter_matrices_ctrl('calculate','10_CD_Ctrl',ei10_C);
%     parameter_matrices_ctrl('calculate','10_CD_Ctrl1',ei10_C1);
    parameter_matrices_ctrl('calculate','10_CC_Ctrl',ei10_A);
    parameter_matrices_ctrl('calculate','10_C_Comb',ei_comb);
    parameter_matrices_ctrl('calculate','10_C_150',ei_150);
%     parameter_matrices('calculate','15_C',ei15_C);
%     parameter_matrices('calculate','10_A',ei10_A);
%     parameter_matrices('calculate','15_A',ei15_A);
    send_email({'samsoon.inayat@uleth.ca'},'Complete - Loading data for Protocol 10 and 15')
catch
    send_email({'samsoon.inayat@uleth.ca'},'Error occurred while loading data')
    lasterror
end
%%
for ii = 1:size(ET15,1)
    ei15(ii) = getData_py(f,ET15(ii,:));
end


%%
ei10 = loadContextsResponses_1(ei10,[1 1],[1 -1 1]);
ei15 = loadContextsResponses_ctrl(ei15,[1 1],[0 0 0]);
ei16 = loadContextsResponses_ctrl(ei16,[1 1],[0 0 0]);

% parameter_matrices('calculate','10_C',ei10_C);
% parameter_matrices('calculate','15_C',ei15_C);
% parameter_matrices('calculate','16',ei16);
disp('All Done!');
%%


ei10 = loadContextsResponses_ctrl(ei10,[0 0],[-1 -1 -1]);
ei10 = loadContextsResponses_ctrl(ei10,[1 1],[0 0 0]);
% ei10 = loadContextsResponses_ctrl(ei10,1,[1 1 1]);

%% for Sam-WS
owr = [1,1]; owrp = [0 0 0];
for ii = 1:length(selRecs)
    if ismember(ii,[8])
        ei10(ii) = loadContextsResponses_ctrl(ei10(ii),owr,owrp);
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
        ei10(ii) = loadContextsResponses_ctrl(ei10(ii),owr,owrp);
    end
end
disp('Done!');

owr = [1,1]; owrp = [1 0 0];
for ii = 1:size(selT,1)
    if ismember(ii,[1:8])
        ei10(ii) = loadContextsResponses_ctrl(ei10(ii),owr,owrp);
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