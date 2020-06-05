% function post_processing_context_definitions(sei)
% behaviorPlot(ei10(8))
sei = ei10;
good = defineContexts_10_first_two_animals(sei(1:3),1);processContextDefinitions(sei(1:3));
good = defineContexts_10_second_two_animals(sei(4:6),1);processContextDefinitions(sei(4:6));
good = defineContexts_10(sei(7:13),1);processContextDefinitions(sei(7:13));

good = defineContexts_10(sei(9),1);processContextDefinitions(sei(9));
% if good
%     processContextDefinitions(sei);
%     disp('Contexts Definitions Done!!!');
% end
