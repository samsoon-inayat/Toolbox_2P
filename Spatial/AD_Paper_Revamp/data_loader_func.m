function figure_place_remapping_AD_trials_all_one

%% Load Data
while 1
    mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size;
    ei_A = evalin('base','ei10_A'); 
    ei_C = evalin('base','ei10_C'); 
    
    selContexts = [1 2 3 4 1 2 3 4 1 2 1 2 3 4];
    rasterNames = {'airD','airD','airD','airD','airIT','airIT','airIT','airIT','beltD','beltD','airT','airT','airT','airT'};
    rasterNamesTxt = {'1-t-D','2-t-D','3-t-D','4-t-D','1-i-T','2-i-T','3-i-T','4-i-T','1-b-D','2-b-D','1-t-T','2-t-T','3-t-T','4-t-T'};
    selContexts = [1 2 3 4 1 2 3 4];
    rasterNames = {'airD','airD','airD','airD','airIT','airIT','airIT','airIT'};
    rasterNamesTxt = {'1-t-D','2-t-D','3-t-D','4-t-D','1-i-T','2-i-T','3-i-T','4-i-T'};
%     selContexts = [1 2 3 4 1 2];
%     rasterNames = {'airD','airD','airD','airD','beltD','beltD'};
%     rasterNamesTxt = {'1-t-D','2-t-D','3-t-D','4-t-D','1-b-D','2-b-D'};
    xlabelsSeq = rasterNamesTxt;

    oA = get_data(ei_A,selContexts,rasterNames);
    oC = get_data(ei_C,selContexts,rasterNames);
    for ii = 1:length(selContexts)
        all_xl{ii} = sprintf('%c%c-%d',rasterNames{ii}(1),rasterNames{ii}(end),selContexts(ii));
    end
    C1_t_D = 1; C2_t_D = 2; C3_t_D = 3; C4_t_D = 4;
    C1_i_T = 5; C2_i_T = 6; C3_i_T = 7; C4_i_T = 8;
    C1_b_D = 9; C2_b_D = 10; 
    C1_t_T = 11; C2_t_T = 12; C3_t_T = 13; C4_t_T = 14;
% %     C1_b_D = 5; C2_b_D = 6; 
    
    break
end
n = 0;
