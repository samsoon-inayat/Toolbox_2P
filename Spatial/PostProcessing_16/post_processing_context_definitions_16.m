% function post_processing_context_definitions(sei)
% behaviorPlot(sei(3))
good = defineContexts_16(sei,1);
if good
    processContextDefinitions(sei);
    disp('Contexts Definitions Done!!!');
end
