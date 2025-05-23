function overall_population_vectors

mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size;
ei_C = evalin('base','ei10_C'); 
ei_A = evalin('base','ei10_A'); 


selContexts = [1 2 3 4];
rasterNames = {'airD','airD','airD','airD'};

ei = ei_C;
Rs = get_rasters_data(ei,selContexts,rasterNames);

mRs = calc_mean_rasters(Rs,1:10);
Rs = find_responsive_rasters(Rs,1:10);
view_population_vector(Rs,mRs,300);

ei = ei_A;
Rs = get_rasters_data(ei,selContexts,rasterNames);

mRs = calc_mean_rasters(Rs,1:10);
Rs = find_responsive_rasters(Rs,1:10);
view_population_vector(Rs,mRs,400);
