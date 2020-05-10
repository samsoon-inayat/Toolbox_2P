function sei = load_contexts_data_15_1(sei)
%% load the contexts into ei variable
sei = loadContextsResponses(sei);
disp('Contexts Loaded!!!');
sei = findSpeedRasters(sei,[1 2 3],0);
%% do gaussian fit on means (for now on distance rasters and for air trials only)
sei = gaussfitOnMeans(sei,'air',0);
sei = gaussfitOnMeansT(sei,'air',0);
sei = gaussfitOnMeansT(sei,'light',0);
disp('Gaussian fitting done!!!');
%% find place cell properties
sei = placeCellProperties(sei,'air',0);
sei = placeCellProperties(sei,'belt',0);
disp('Place cell properties done!!!');
