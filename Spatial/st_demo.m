function st_demo
ei = evalin('base','ei{1}');

% Demo program
% STeP is applied to simulated data
%
% Copyright (C) 2016, Yusuke Takeda, ATR, takeda@atr.jp

% Set parameters for this simulation test
% 
% clear all
% close all
% 
% signals = ei.tP.signals;
% 
% data = signals(ei.areCells,:)';

for ii = 1:length(ei.deconv.spSigAll)
    allData(ii,:) = ei.deconv.spSigAll{ii};
end


b = ei.b;
% onset = b.photo_sensor_f(1:(length(b.photo_sensor_f)-1));
% offset = b.photo_sensor_f(2:end);
onset = b.air_puff_r(b.sTrials);
offset = b.air_puff_f(b.sTrials)+100;
trialNumbers = 1:length(onset);
data = [];

for iii = 1:7%length(trialNumbers)
    ii = trialNumbers(iii);
    st = onset(ii);
    se = offset(ii);
%     se = offset(ii);
    frames = (find(b.frames_f >= st & b.frames_f <= se)); % find frame numbers
    if iii == 1
        onsetM(iii) = 1;
        offsetM(iii) = length(frames);
        firstFrame = frames(1);
    else
        onsetM(iii) = size(data,2)+1;
        offsetM(iii) = size(data,2)+length(frames);
    end
    data = [data allData(:,frames)];
    st = offset(ii);
    se = onset(ii+1);
%     se = offset(ii);
    frames = (find(b.frames_f >= st & b.frames_f <= se)); % find frame numbers
    data = [data allData(:,frames)];
end


data = data';

data_length=size(data,1);% Length of simulated data
N=60;% Length of spatiotemporal pattern
Npattern=10;% Number of spatiotemporal patterns

%% Estimate spatiotemporal patterns and their onsets using STeP

fileName = 'sResults.mat';%makeName(sprintf('stResults.mat',ii,N,Npattern),ei.folders.thispFolder);
if exist(fileName,'file')
    load(fileName);
else
%     Apply STeP to simulated data
    [eonset,epattern]=st_STeP(data(1:data_length,:),N,Npattern);
    save(fileName,'eonset','epattern');
end


ii = ei.db.selectedPlane;
allFrames = ii:(ei.zSteps+1):(ei.totalFrames*(ei.zSteps+1));
% Show estimation result
eonset_timeseries=st_make_onset_timeseries(eonset,data_length);
for pt=1:Npattern
    thisPattern = squeeze(epattern(:,pt,:))';
    thisOnsets = eonset_timeseries(:,pt);
    maxThisPattern = max(thisPattern,[],2);
    [~,sinds] = sort(maxThisPattern);
    sinds = flipud(sinds);
    centroids = getCentroids(ei,sinds(1:10));
    sThisPattern = thisPattern(sinds,:);
    figure(101);clf;
    imagesc(sThisPattern);
    ops = find(thisOnsets);
    dops = find(diff(ops)<10);
    ops(dops) = [];
    for eo = 1:length(ops)
        tO = ops(eo)+firstFrame;
        thisFrames = tO:(tO+9);
        frameNums = allFrames(thisFrames);
        frames = getFramesFromRaw(ei,frameNums);
        rFrames = register_movie(frames,ei.ops1{1},ei.ops1{1}.DS(thisFrames,:));
        for ff = 1:size(frames,3)
            figure(102);clf;
            imagesc(rFrames(:,:,ff));
            for mm = 1:10
                text(centroids(mm,1),centroids(mm,2),num2str(sinds(mm)));
            end
            titleText = sprintf('%d',tO);
            title(titleText);
            pause(0.1);
        end
        pause;
    end
    n = 0;
end

return;

figure
for pt=1:Npattern
    subplot(Npattern,5,5*(pt-1)+1)
    imagesc(squeeze(epattern(:,pt,:))')
    colorbar
    if pt==1
        title('Estimated spatiotemporal pattern')
    end
    if pt==Npattern
        xlabel('Time')
    end
    ylabel('Channel')
    subplot(Npattern,5,5*(pt-1)+2:5*pt)
    plot(eonset_timeseries(:,pt));hold on;
    plot(onsetM,ones(size(onsetM))*0.5,'ro');
    plot(offsetM,ones(size(offsetM))*0.25,'mo');
    if pt==1
        title('Estimated onset timeseries')
    end
    if pt==Npattern
        xlabel('Time')
    end
end



