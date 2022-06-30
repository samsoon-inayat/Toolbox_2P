function figure_place_remapping_AD_trials_all_one

%% Load Data
while 1
    mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size;
    ei_A = evalin('base','ei10_A'); 
    ei_C = evalin('base','ei10_C'); 
    
    selContexts = [1 2 3 4 1 2 3 4];
    rasterNames = {'airD','airD','airD','airD','airIT','airIT','airIT','airIT'};
    rasterNamesTxt = {'1-D','2-D','3-D','4-D','1-T','2-T','3-T','4-T'};
    xlabelsSeq = rasterNamesTxt;

    oA = get_data(ei_A,selContexts,rasterNames);
    oC = get_data(ei_C,selContexts,rasterNames);
    for ii = 1:length(selContexts)
        all_xl{ii} = sprintf('%c%c-%d',rasterNames{ii}(1),rasterNames{ii}(end),selContexts(ii));
    end
    C1_t_D = 1; C2_t_D = 2; C3_t_D = 3; C4_t_D = 4;
    C1_i_T = 5; C2_i_T = 6; C3_i_T = 7; C4_i_T = 8;
    
    break
end
n = 0;

