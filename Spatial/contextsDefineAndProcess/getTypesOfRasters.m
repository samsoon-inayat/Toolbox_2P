function typesOfRasters = getTypesOfRasters(allContexts,thisStimMarker)

ind = strcmp(allContexts.typesOfMarkers,thisStimMarker);

typesOfRasters = allContexts.typesOfRasters{ind};

if ischar(typesOfRasters)
    typesOfRasters = {typesOfRasters};
end