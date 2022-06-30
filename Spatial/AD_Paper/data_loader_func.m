function figure_place_remapping_AD_trials_all_one

%% Load Data
while 1
    mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size;
    ei_A = evalin('base','ei10_A'); 
    ei_C = evalin('base','ei10_C'); 
    
    selContexts = [1 2 3 4];
    rasterNames = {'airD','airD','airD','airD'};
    rasterNamesTxt = {'C1','C2','C3','C4'};
    xlabelsSeq = rasterNamesTxt;

    oA = get_data(ei_A,selContexts,rasterNames);
    oC = get_data(ei_C,selContexts,rasterNames);
    for ii = 1:length(selContexts)
        all_xl{ii} = sprintf('%c%c-%d',rasterNames{ii}(1),rasterNames{ii}(end),selContexts(ii));
    end
    C1 = 1; C2 = 2; C3 = 3; C4 = 4;
    
    break
end
n = 0;

%%
props1 = get_props_Rs(o.Rs,50);