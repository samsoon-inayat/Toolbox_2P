function conjunctive_representation

%% Load Data
tic
mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size;
ei = evalin('base','ei');
udata = evalin('base','udata1');
% udata1 = evalin('base','udata1');
% udataT = evalin('base','udataT');
% udataD = evalin('base','udataD');

outT = get_the_binned_data(udata,'time',0.3);
outD = get_the_binned_data(udata,'dist',3);

filename = fullfile(mData.pd_folder,'FR_tuning.mat');
if exist(filename,'file')
    load(filename);
else
    % Get the metrics and run stats
    variable_combs = {'FR_time','FR_dist','FR_speed'};
    trialsOrConcat = 'concatenate'; %'trials' or "concatenate"
    nshuffles = 1000;
    clc
    for vn = 1:length(variable_combs)
        vn
        met_valsT{vn} = get_metrics_FR(outT,variable_combs{vn},trialsOrConcat,nshuffles);
        met_valsD{vn} = get_metrics_FR(outD,variable_combs{vn},trialsOrConcat,nshuffles);
    end
    
    save(filename,"met_valsT","met_valsD");
end

% delete second plane data
for ii = 1:5
    if length(ei{ii}.plane) > 1
        ei{ii}.plane(2) = [];
    end
    disp(length(ei{ii}.plane))
end

% [udata,udata1] = get_unlv_analysis_data(ei);
% udataT = get_unlv_analysis_data_time_bin(udata1,0.3);
% udataD = get_unlv_analysis_data_dist_bin(udata1,3);


selContexts = [1 4 6 2 7 3 3 4 4 5 5 3 3 4 4 5 5 0 0 2 2 7 7 2 7 2 7 3 4 5 3 4 5 3 4 5 3 4 5];
rasterNames = {'light22T','light22T','light22T','air55T','air55T','airD','airID','airD','airID','airD','airID','airT','airIT','airT','airIT','airT','airIT','motionOnsets','motionOffsets','airRT','airIRT','airRT','airIRT',...
    'airOnsets22T','airOnsets22T','airOffsets22T',...
                'airOffsets22T','airOnsets22T','airOnsets22T','airOnsets22T',...
                'airOffsets22T','airOffsets22T','airOffsets22T','beltD','beltD','beltD','beltT','beltT','beltT'};
rasterNamesTxt = {'Lb','ArL-L','Lb*','Ab','Ab*','C3-AOn-D','C3-AOff-D','C4-AOn-D','C4-AOff-D','C5-AOn-D','C5-AOff-D','C3-AOn','C3-AOff','C4-AOn','C4-AOff','C5-AOn','C5-AOff','MOn','MOff','C2-AOn','C2-AOff','C7-AOn','C7-AOff',...
    '2-AOn','7-AOn','2-AOff','7-AOff',...
    '3-AOn','4-AOn','5-AOn','3-AOff','4-AOff','5-AOff','3-B-D','4-B-D','5-B-D','3-B','4-B','5-B'};
% rasterNamesTxt = {'Lb','ArL-L','Lb*','Ab','Ab*','Ar-t-D','Ar-i-D','ArL-t-D','ArL-i-D','Ar*-t-D','Ar*-i-D','Ar-t','Ar-i','ArL-t','ArL-i','Ar*-t','Ar*-i','MOn','MOff','Ab-t','Ab-i','Ab*-t','Ab*-i',...
%     '2-AOn','7-AOn','2-AOff','7-AOff',...
%     '3-AOn','4-AOn','5-AOn','3-AOff','4-AOff','5-AOff','3-B-D','4-B-D','5-B-D','3-B','4-B','5-B'};
xlabelsSeq = {'Lb','ArL-L','Lb*','Ab','Ab*','Ar-t','ArL-t','Ar*-t','Ar-i','ArL-i','Ar*-i'};

o = get_data(ei,selContexts,rasterNames);
for ii = 1:length(selContexts)
    all_xl{ii} = sprintf('%c%c-%d',rasterNames{ii}(1),rasterNames{ii}(end),selContexts(ii));
end
Lb_T = 1; ArL_L_T = 2; Lbs_T = 3; Ab_T = 4; Abs_T = 5; Ab_t_T = 20; Abs_t_T = 22;  Ab_i_T = 21; Abs_i_T = 23;
Ar_t_D = 6; Ar_t_T = 12; ArL_t_D = 8; ArL_t_T = 14; Ars_t_D = 10; Ars_t_T = 16;
Ar_i_D = 7; Ar_i_T = 13; ArL_i_D = 9; ArL_i_T = 15; Ars_i_D = 11; Ars_i_T = 17;
MOn_T = 18; MOff_T = 19;
Ab_On = 24; Abs_On = 25; Ab_Off = 26; Abs_Off = 27; 
Ar_On = 28; ArL_On = 29; Ars_On = 30; Ar_Off = 31; ArL_Off = 32; Ars_Off = 33; Ar_B_D = 34; ArL_B_D = 35; Ars_B_D = 36;
Ar_B_T = 37; ArL_B_T = 38; Ars_B_T = 39;

    [speedRs,resp_speed,speed_percent,resp_speedAcc] = load_speed_response(ei);
    % [speedRsZ] = load_speed_responseZ(ei);
%     all_xl{ii+1} = 'sp';
%     resp = [o.resp.vals resp_speed];

M = [18 19];

%     dzMI = prop_op(o.props.zMI(:,[Ar_t_D Ar_i_D]),o.props.zMI(:,[Ar_t_T Ar_i_T]),0.1);

toc
%
for rr = 1:size(o.Rs,1)
    for cc = 1:size(o.Rs,2)
        t_resp = o.Rs{rr,cc}.resp;
        if length(t_resp.trial_scores) < length(t_resp.vals)
            num_cells = length(t_resp.trial_scores);
%         if isfield(t_resp,'trial_scores')
%             disp(length(t_resp.trial_scores));
            fields = fieldnames(t_resp);
            for ii = 1:length(fields)
                cmdTxt = sprintf('sz = (size(t_resp.%s));',fields{ii}); 
                eval(cmdTxt);
                if sum(sz>num_cells) > 0
                    if sz(1) > num_cells
                        cmdTxt = sprintf('t_resp.%s = t_resp.%s(1:num_cells,:);',fields{ii},fields{ii}); 
                    else
                        cmdTxt = sprintf('t_resp.%s = t_resp.%s(:,1:num_cells);',fields{ii},fields{ii}); 
                    end
                    eval(cmdTxt);
                end
            end
            o.Rs{rr,cc}.resp = t_resp;
        end
    end
end

for ii = 1:length(ei)
    centroids{ii} = ei{1}.plane{1}.tP.centroids;
    dists_centroids{ii} = ei{1}.plane{1}.tP.dists_centroids;
end

cellposeCells = evalin('base','cellposeCells');

o = get_trials_MI(o);
n = 0;
n = 0;
n = 0;
n = 0;
n = 0;
%%
% %% Load Data
% tic
% mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size;
% ei = evalin('base','ei'); 
% 
% selContexts = [3 3 4 4 5 5];
% rasterNames = {'airT','airIT','airT','airIT','airT','airIT'};
% rasterNamesTxt = {'Lb','ArL-L','Lb*','Ab','Ab*','Ar-t-D','Ar-i-D','ArL-t-D','ArL-i-D','Ar*-t-D','Ar*-i-D','Ar-t','Ar-i','ArL-t','ArL-i','Ar*-t','Ar*-i','MOn','MOff','Ab-t','Ab-i','Ab*-t','Ab*-i',...
%     '2-AOn','7-AOn','2-AOff','7-AOff',...
%     '3-AOn','4-AOn','5-AOn','3-AOff','4-AOff','5-AOff','3-B-D','4-B-D','5-B-D','3-B','4-B','5-B'};
% xlabelsSeq = {'Lb','ArL-L','Lb*','Ab','Ab*','Ar-t','ArL-t','Ar*-t','Ar-i','ArL-i','Ar*-i'};
% 
% o = get_data(ei,selContexts,rasterNames);
% for ii = 1:length(selContexts)
%     all_xl{ii} = sprintf('%c%c-%d',rasterNames{ii}(1),rasterNames{ii}(end),selContexts(ii));
% end
% Lb_T = 1; ArL_L_T = 2; Lbs_T = 3; Ab_T = 4; Abs_T = 5; Ab_t_T = 20; Abs_t_T = 22;  Ab_i_T = 21; Abs_i_T = 23;
% Ar_t_D = 6; Ar_t_T = 12; ArL_t_D = 8; ArL_t_T = 14; Ars_t_D = 10; Ars_t_T = 16;
% Ar_i_D = 7; Ar_i_T = 13; ArL_i_D = 9; ArL_i_T = 15; Ars_i_D = 11; Ars_i_T = 17;
% MOn_T = 18; MOff_T = 19;
% Ab_On = 24; Abs_On = 25; Ab_Off = 26; Abs_Off = 27; 
% Ar_On = 28; ArL_On = 29; Ars_On = 30; Ar_Off = 31; ArL_Off = 32; Ars_Off = 33; Ar_B_D = 34; ArL_B_D = 35; Ars_B_D = 36;
% Ar_B_T = 37; ArL_B_T = 38; Ars_B_T = 39;
% 
% %     [speedRs,resp_speed,speed_percent,resp_speedAcc] = load_speed_response(ei);
% %     all_xl{ii+1} = 'sp';
% %     resp = [o.resp.vals resp_speed];
% 
% M = [18 19];
% 
% %     dzMI = prop_op(o.props.zMI(:,[Ar_t_D Ar_i_D]),o.props.zMI(:,[Ar_t_T Ar_i_T]),0.1);
% 
% toc
% n = 0;

