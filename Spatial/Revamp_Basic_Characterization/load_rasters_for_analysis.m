function load_rasters_for_analysis

%% Load Data
while 1
    mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size;
    ei = evalin('base','ei'); 
    
    selContexts = [3 4 5];
    rasterNames = {'air77T','air77T','air77T'};
    o = get_data(ei,selContexts,rasterNames);
    break
end
n = 0;
%%