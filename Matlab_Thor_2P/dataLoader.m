%%
add_to_path
%% load data
clear all
clc

owr = 0;
ei{1} = getAllData('173062','2018-07-13',1,''); ei{2} = getAllData('173062','2018-07-13',2,'');
ei{3} = getAllData('173511','2018-07-13',1,''); ei{4} = getAllData('173511','2018-07-13',2,'');
ei{5} = getAllData('173198','2018-07-13',1,''); ei{6} = getAllData('173198','2018-07-13',2,'');

% ei{1} = getAllData('173062','2018-06-22',1,''); ei{2} = getAllData('173062','2018-06-22',2,'');
% ei{1} = getAllData('173062','2018-07-12',1,''); ei{2} = getAllData('173062','2018-07-12',2,'');
% ei{3} = getAllData('173511','2018-07-12',1,''); ei{4} = getAllData('173511','2018-07-12',2,'');

ei{7} = getAllData('173062','2018-07-24',1,''); ei{8} = getAllData('173062','2018-07-24',2,'');

for ii = [5:8]
    ei(ii) = loadContextsResponses(ei(ii));
end
