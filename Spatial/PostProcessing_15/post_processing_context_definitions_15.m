% function post_processing_context_definitions(sei)
% behaviorPlot(sei(1))
good = defineContexts_protocol15(sei,1);
if good
    processContextDefinitions(sei);
    disp('Contexts Definitions Done!!!');
end
