

sei = [];
sei = ei(1);%(4:5);%([4 5]);
disp('Done');


% for ii = 1:length(sei)
%     for jj = 1:length(sei{ii}.plane)
%         disp([ii length(sei{ii}.plane{jj}.b.frames_f) length(sei{ii}.plane{jj}.tP.deconv.spSigAll{1})])
%     end
% end 

% for ii = 1:7
%     disp(ei{1}.plane{1}.contexts(ii).trials);
% end

%%
% After selecting ei to work with, if context definitions need to be
% processed run
% post_processing_context_definitions_15
% sei = process_contexts_15(sei)
% processContexts(ei15,1);
%%
ei15 = findSpeedRasters(ei15,[3 4 5],0);
ei15 = load_contexts_data_15(ei15);

%% get data of responses only that is only the variable rasters from ei
[data15,mData15] = getRasterData(ei15,[3 4 5],'air');
[data15b,mData15b] = getRasterData(ei15,[3 4 5],'belt');
[data15AOn,mData15AOn] = getRasterData(ei15,[3 4 5],'airOnsets27');
[data15AOff,mData15AOff] = getRasterData(ei15,[3 4 5],'airOffsets27');
[data15L,mData15L] = getRasterData(ei15,[1 6],'light');
[data15AOnR,mData15AOffR] = getRasterData(ei15,[2 7],'airOnsets11');
[data15AOffR,mData15AOffR] = getRasterData(ei15,[2 7],'airOffsets11');

% mData.colors = {[0 0 0],[0.1 0.7 0.l3],'r','b','m','c','g','y'};
mData15.colors = {'k','b','r','m','c','g','y'};
% mData.colors = getColors(10,{'w','g'});
mData15.axes_font_size = 12;
mData15.sigColor = [0.54 0.27 0.06];
mData15.pdf_folder = fullfile(pwd,'PDFs');
disp('data extracted');