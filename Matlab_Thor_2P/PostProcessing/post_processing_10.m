

sei = [];
sei = eip;%(4:5);%([4 5]);
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
% sei = process_contexts(sei)
%%
sei = findSpeedRasters(sei,[1 2],1);
%%
eip = findSpeedRasters(eip,[1 2],1);
eip = load_contexts_data_10(eip);
%% get data of responses only that is only the variable rasters from ei
data = []; datab = [];
[data,mData] = getRasterData(eip,[1 2],'air');
[datab,mDatab] = getRasterData(eip,[3 4],'belt');
for ii = 1:length(eip)
    for pp = 1:length(data{ii})
%         disp([ii pp])
        data{ii}{pp}{3} = datab{ii}{pp}{3-2};
        data{ii}{pp}{4} = datab{ii}{pp}{4-2};
    end
end
mData.colors = {[0 0 0],[0.1 0.7 0.3],'r','b','m','c','g','y'};
% mData.colors = getColors(10,{'w','g'});
mData.axes_font_size = 12;
mData.sigColor = [0.54 0.27 0.06];
mData.pdf_folder = fullfile(pwd,'PDFs');
disp('data extracted');