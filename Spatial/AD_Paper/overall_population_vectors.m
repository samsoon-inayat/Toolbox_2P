function overall_population_vectors

mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size;
ei_C = evalin('base','ei10_C2'); 
ei_A = evalin('base','ei10_A2'); 


selContexts = [1 2 3 4];
rasterNames = {'airD','airD','airD','airD'};

ei = ei_C;
Rs = get_rasters_data(ei,selContexts,rasterNames);

mRs = calc_mean_rasters(Rs,1:10);
Rs = find_responsive_rasters(Rs,1:10);
view_population_vector(Rs,mRs,300);
% view_population_vector_corr(Rs,mRs,400);
% [CR,aCR] = find_population_vector_corr(Rs,mRs,650);

ei = ei_A;
Rs = get_rasters_data(ei,selContexts,rasterNames);

mRs = calc_mean_rasters(Rs,1:10);
Rs = find_responsive_rasters(Rs,1:10);
view_population_vector(Rs,mRs,500);
% view_population_vector_corr(Rs,mRs,600);
% [CR,aCR] = find_population_vector_corr(Rs,mRs,1);
% [CR,aCR] = find_population_vector_corr(Rs([1 3 4 5],:),mRs([1 3 4 5],:),700);


