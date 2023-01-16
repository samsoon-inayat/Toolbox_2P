function conjunctive_representation

% %% Load Data
% tic
% mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size;
% ei = evalin('base','ei'); 
% 
% selContexts = [1 4 6 2 7 3 3 4 4 5 5 3 3 4 4 5 5 0 0 2 2 7 7 2 7 2 7 3 4 5 3 4 5 3 4 5 3 4 5];
% rasterNames = {'light22T','light22T','light22T','air55T','air55T','airD','airID','airD','airID','airD','airID','airT','airIT','airT','airIT','airT','airIT','motionOnsets','motionOffsets','airRT','airIRT','airRT','airIRT',...
%     'airOnsets22T','airOnsets22T','airOffsets22T',...
%                 'airOffsets22T','airOnsets22T','airOnsets22T','airOnsets22T',...
%                 'airOffsets22T','airOffsets22T','airOffsets22T','beltD','beltD','beltD','beltT','beltT','beltT'};
% rasterNamesTxt = {'Lb-T','ArL-L-T','Lb*-T','Ab-T','Ab*-T','Ar-t-D','Ar-i-D','ArL-t-D','ArL-i-D','Ar*-t-D','Ar*-i-D','Ar-t-T','Ar-i-T','ArL-t-T','ArL-i-T','Ar*-t-T','Ar*-i-T','MOn-T','MOff-T','Ab-t-T','Ab-i-T','Ab*-t-T','Ab*-i-T',...
%     '2-AOn','7-AOn','2-AOff','7-AOff',...
%     '3-AOn','4-AOn','5-AOn','3-AOff','4-AOff','5-AOff','3-B-D','4-B-D','5-B-D','3-B-T','4-B-T','5-B-T'};
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

%% Load Data
tic
mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size;
ei = evalin('base','ei'); 

selContexts = [3 3 4 4 5 5];
rasterNames = {'airT','airIT','airT','airIT','airT','airIT'};
rasterNamesTxt = {'Lb-T','ArL-L-T','Lb*-T','Ab-T','Ab*-T','Ar-t-D','Ar-i-D','ArL-t-D','ArL-i-D','Ar*-t-D','Ar*-i-D','Ar-t-T','Ar-i-T','ArL-t-T','ArL-i-T','Ar*-t-T','Ar*-i-T','MOn-T','MOff-T','Ab-t-T','Ab-i-T','Ab*-t-T','Ab*-i-T',...
    '2-AOn','7-AOn','2-AOff','7-AOff',...
    '3-AOn','4-AOn','5-AOn','3-AOff','4-AOff','5-AOff','3-B-D','4-B-D','5-B-D','3-B-T','4-B-T','5-B-T'};
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

%     [speedRs,resp_speed,speed_percent,resp_speedAcc] = load_speed_response(ei);
%     all_xl{ii+1} = 'sp';
%     resp = [o.resp.vals resp_speed];

M = [18 19];

%     dzMI = prop_op(o.props.zMI(:,[Ar_t_D Ar_i_D]),o.props.zMI(:,[Ar_t_T Ar_i_T]),0.1);

toc
n = 0;

