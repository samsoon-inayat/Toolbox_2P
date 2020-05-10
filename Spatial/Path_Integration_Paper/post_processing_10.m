

sei = [];
sei = ei10;%(4:5);%([4 5]);
disp('Done');


% for ii = 1:length(sei)
%     for jj = 1:length(sei{ii}.plane)
%         disp([ii length(sei{ii}.plane{jj}.b.frames_f) length(sei{ii}.plane{jj}.tP.deconv.spSigAll{1})])
%     end
% end 

%%
% After selecting ei to work with, if context definitions need to be
% processed run

% post_processing_context_definitions
% 
% % behaviorPlot(ei10(3))
% good = defineContexts_10_first_two_animals(ei10(1:2),1);processContextDefinitions(ei10(1:2));
% good = defineContexts_10_second_two_animals(ei10(3:4),1);processContextDefinitions(ei10(3:4));
% good = defineContexts_10(ei10(5:8),1);processContextDefinitions(ei10(5:8));
% find_trial_distance_time_distributions(ei10,0)

% ei10 = loadContextsResponses(ei10,0,[0 0 0]);
% ei10 = correctingErrors(ei10,0,[0 1 1]);
ei10 = getData_py(f,T10.T([7 11 13 14],:));
ei10 = loadContextsResponses(ei10(1:2),-1,[-1 -1 0]);


owr = -1; owrp = [-1 1 -1];
for ii = 1:size(selT,1)
    if ismember(ii,[1:8])
        ei10(ii) = loadContextsResponses(ei10(ii),owr,owrp);
    end
end

%%
mData.colors = {[0 0 0],[0.1 0.7 0.3],'r','b','m','c','g','y'};
% mData.colors = getColors(10,{'w','g'});
mData.axes_font_size = 12;
mData.sigColor = [0.54 0.27 0.06];
mData.pdf_folder = fullfile(pwd,'PDFs');
disp('data extracted');