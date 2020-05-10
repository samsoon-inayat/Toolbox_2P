

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
%%
% ei = findSpeedRasters(ei,[3 4 5],1);

ei = load_contexts_data_15(ei);

%% get data of responses only that is only the variable rasters from ei
data = []; dataL = []; dataA =[];
[data,mData] = getRasterData(ei,[3 4 5],'air');
[dataL,mDataL] = getRasterDataT(ei,[1 6],'light');
[dataA,mDataA] = getRasterDataT(ei,[2 7],'air');

% mData.colors = {[0 0 0],[0.1 0.7 0.3],'r','b','m','c','g','y'};
mData.colors = {'k','b','r','m','c','g','y'};
% mData.colors = getColors(10,{'w','g'});
mData.axes_font_size = 12;
mData.sigColor = [0.54 0.27 0.06];
mData.pdf_folder = fullfile(pwd,'PDFs');
disp('data extracted');