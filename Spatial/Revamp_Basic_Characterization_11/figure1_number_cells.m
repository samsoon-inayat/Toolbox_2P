function figure1_number_PCs_threshold_comparison

%% Load Data
while 1
    mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size;
    ei = evalin('base','ei'); 
    
    selContexts = [3 3 4 4 5 5]; rasterNames = {'airD','airIT','airD','airIT','airD','airIT'};
    oDT = get_data(ei,selContexts,rasterNames);

    break
end
n = 0;
%%
if 1
    RsC = oDT.Rs;
    for rr = 1:size(RsC,1)
        numset1(rr) = length(RsC{rr,1}.iscell);
    end
    return;
end

