function sample_video

mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size;
ei = evalin('base','ei'); 

selContexts = [1];
rasterNames = {'airD'};
% Rs = get_rasters_data(ei_11_15,selContexts,rasterNames);
Rs = get_rasters_data(ei,selContexts,rasterNames);
Rs = find_responsive_rasters(Rs,1:10);
mR = calc_mean_rasters(Rs,1:10);
% [CRc,aCRc,mRR] = find_population_vector_corr(Rs,mR,1);
% [resp_fractionC,resp_valsC,OIC,mean_OIC,resp_ORC,resp_OR_fractionC,resp_ANDC,resp_AND_fractionC] = get_responsive_fraction(Rs);
n = 0;

%%

vo = evalin('base','vo');
trigger_frame = evalin('base','trigger_frame');
%%
currentTime = vo.CurrentTime;
%%
time_of_trigger_frame = trigger_frame/vo.FrameRate;
vo.CurrentTime = time_of_trigger_frame;

%%
b = ei{1}.b;
b.frames_f = ei{1}.plane{1}.b.frames_f;
bts = b.ts;
vfr = trigger_frame:floor((vo.Duration*60));
vts = (vfr-trigger_frame)/vo.FrameRate;
n = 0;
%%

trialStart = 20;
trialEnd = 20;

tStart = bts(b.air_puff_r(trialStart)-round(1e6 * 1/b.si));
tEnd = bts(b.air_puff_f(trialEnd)+round(1e6 * 1/b.si));
vfrStart = get_video_frame_and_time(vfr,vts,tStart);
vfrEnd = get_video_frame_and_time(vfr,vts,tEnd);


trials = trialStart:trialEnd
for ii = 1:length(trials)
    tStart1 = bts(b.air_puff_r(trials(ii)));
    tEnd1 = bts(b.air_puff_f(trials(ii)));
    vfrStartT(ii) = get_video_frame_and_time(vfr,vts,tStart1);
    vfrEndT(ii) = get_video_frame_and_time(vfr,vts,tEnd1);
end

%% get calcium imaging video
indS = find(b.ts >= tStart,1,'first'); indE = find(b.ts >= tEnd,1,'first');
frameNumbers = find(b.frames_f >= indS & b.frames_f <= indE);
rows = 512; cols = 512;
rawFile = fullfile(ei{1}.plane{1}.s2p_folder,'data.bin');
rawFile = ei{1}.thorExp.rawFile;
frames = getFramesFromRaw_simple (rawFile,rows,cols,frameNumbers);
sT = evalin('base','T_C');
average_image = max(frames,[],3);
% navgimg = normalize_image(average_image);

%% get behavior video
frameNums = vfrStart:vfrEnd;
vFrames = get_video_frames(vo,frameNums);

%%
m = min([min(average_image(:)) min(frames(:))])+100;
M = 0.5*max([max(average_image(:)) max(frames(:))]);
figure(1000);clf;
subplot 211;
hi = imagesc(vFrames{1});
axis off;
axis equal;
subplot 212;
ind = 1;
cFrame = imfuse((average_image),frames(:,:,ind),'blend');
cFrame = imfuse((average_image),cFrame,'blend');
hic = imagesc((cFrame)); ind = ind + 1;
axis off;
axis equal;
colormap gray;
title(frameNums(1));

for ii = 2:length(frameNums)
    set(hi,'CData',vFrames{ii});
    if ismember(frameNums(ii), vfrStartT)
        hts = text(260,260,'Air','color','r');
    end
    if ismember(frameNums(ii),vfrEndT)
        delete(hts);
        break;
    end
    title(frameNums(ii));
    if mod(ii,2) == 0
        cFrame = imfuse((average_image),frames(:,:,ind),'blend');
        cFrame = imfuse((average_image),cFrame,'blend');
        set(hic,'CData',(cFrame));
%         set(hic,'clims',[1000 2500]);
%         xlim([257 512]);
%         ylim([257 512]);
        ind = ind+1;
        colormap gray;
    end
    pause(0.01);
end


function [f,t] = get_video_frame_and_time(vfr,vts,tStart)
ind = find(vts-tStart > 0,1,'first');
f = vfr(ind);
t = vts(ind);

function vFrames = get_video_frames(vo,frameNums)
% frameNums = vfrStart:vfrEnd;
vFrames = [];
for ii = 1:length(frameNums)
    vFrames{ii} = read(vo,frameNums(ii));
end