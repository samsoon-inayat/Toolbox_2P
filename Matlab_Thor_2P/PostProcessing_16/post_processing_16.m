ei

sei = [];
sei = ei;%(1);%(4:5);%([4 5]);
disp('Done');


% for ii = 1:length(sei)
%     for jj = 1:length(sei{ii}.plane)
%         disp([ii length(sei{ii}.plane{jj}.b.frames_f) length(sei{ii}.plane{jj}.tP.deconv.spSigAll{1})])
%     end
% end 

%%
% After selecting ei to work with, if context definitions need to be
% processed run
% post_processing_context_definitions_16
% sei = process_contexts(sei)
%%
ei = load_contexts_data_16(ei);


% for ii = 1:10
%     disp(ei{1}.plane{1}.contexts(ii).trials);
% end
%% get data of responses only that is only the variable rasters from ei
data = []; dataT = []; dataL = []; dataA = [];
[data,mData] = getRasterData(ei,[1 2 3 4 5 6 7],'belt');
[dataT,mDataT] = getRasterDataT(ei,[8],'tone');
[dataL,mDataL] = getRasterDataT(ei,[9],'light');
[dataA,mDataA] = getRasterDataT(ei,[10],'air');
% % [datab,mDatab] = getRasterData(ei,[3 4],'belt');
% for ii = 1:length(ei)
%     for pp = 1:length(data{ii})
% %         disp([ii pp])
%         data{ii}{pp}{3} = datab{ii}{pp}{3-2};
%         data{ii}{pp}{4} = datab{ii}{pp}{4-2};
%     end
% end
mData.colors = {[0 0 0],[0.1 0.7 0.3],'r','b','m','c','g','y'};
% mData.colors = getColors(10,{'w','r','g','m','y'});
mData.axes_font_size = 12;
mData.sigColor = [0.54 0.27 0.06];
mData.pdf_folder = fullfile(pwd,'PDFs');
disp('data extracted');