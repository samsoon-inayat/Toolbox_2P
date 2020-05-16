% function post_processing_context_definitions(sei)
% behaviorPlot(sei(9))
good = defineContexts_10_first_two_animals(sei(1:2),1);processContextDefinitions(sei(1:2));
good = defineContexts_10_second_two_animals(sei(3:4),1);processContextDefinitions(sei(3:4));
good = defineContexts_10(sei(5:8),1);processContextDefinitions(sei(5:8));

good = defineContexts_10(sei(9),1);processContextDefinitions(sei(9));
% if good
%     processContextDefinitions(sei);
%     disp('Contexts Definitions Done!!!');
% end
