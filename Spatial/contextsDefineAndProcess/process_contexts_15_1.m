function sei = process_contexts_15(sei)

processContexts(sei,0);
disp('Context processing done!!!');
%% load the contexts into ei variable
sei = loadContextsResponses_15(sei);
disp('Contexts Loaded!!!');
%% do gaussian fit on means (for now on distance rasters and for air trials only)
% sei = gaussfitOnMeans(sei,'belt',0);
sei = gaussfitOnMeansT(sei,'air',0);
% sei = gaussfitOnMeansT(sei,'light',0);
disp('Gaussian fitting done!!!');
%% find place cell properties
sei = placeCellProperties(sei,'air',0);
% sei = placeCellProperties(sei,'belt',0);
disp('Place cell properties done!!!');
