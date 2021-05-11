function allContexts = contextDefinitions

allContexts.names = {
    'Light - Brake'; % stop the animal with the brake and apply flashes of light
    'Air - Brake'; % stop the animal with the brake and apply air puff
    'Air - No Cues - No Brake';
    'Air - Light - No Brake';
    'Air - Cues - No Brake';
    'Air - Cues - Brake';
    'Air - No Cues - Brake';
    'Air - Tone - Brake';
    'Air - Light - Brake';
    'Tone - Brake'; % stop the animal with the brake and apply tone
    'Spontaneous - No Brake';
    'Spontaneous - Brake';
    };


allContexts.ids = (1:length(allContexts.names))';

% display the following table to fill the variable contextIDs
% such that you are letting the program know the sequence of contexts the 
% mouse was exposed to
table(allContexts.ids,allContexts.names);

allContexts.typesOfMarkers = {'air';'belt';'motionOnsetsOffsets';'airOnsets22';'airOffsets22';'motionI';'motionOnsets22';'motionOffsets22';...
    'airI';'motionOffsetAirOnset';'light';'tone';'airOnsets27';'airOffsets27';'airOnsets11';'airOffsets11';'airOnsets010';'airOffsets010';...
    'airOnsets01';'airOffsets01';'light11';'tone11';'light22';'air44';'tone22';'light0p30p3';'tone0p30p3';'air33';'air77';'air55'};

allContexts.typesOfRasters = {{'dist';'time'};{'dist';'time'};{'dist';'time'};'time';'time';'time';'time';'time';...
    {'dist';'time'};'time';'time';'time';'time';'time';'time';'time';'time';'time';...
    'time';'time';'time';'time';'time';'time';'time';'time';'time';'time';'time';'time'};

allContexts.rasterFunctions = {'find_MI';'getFits_myGaussFit';'fractalDim';'find_MI_1';'getFits_myGaussFit_1';'fractalDim_1';'find_MI_2';'getFits_myGaussFit_2';'fractalDim_2'};
allContexts.varNames = {'info_metrics';'gauss_fit_on_mean';'fractal_dim';'info_metrics_1';'gauss_fit_on_mean_1';'fractal_dim_1';'info_metrics_2';'gauss_fit_on_mean_2';'fractal_dim_2'};



