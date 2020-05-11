function getSpeedResponsePlot(ei)

if ~exist('ei','var')
    aei = evalin('base','ei');
    ei = aei([1 2]);
end

cei = combineRecordings(ei);

b = ei{1}.b;
tei = ei{1};
[a_ddata,a_tdata] = getDataContexts(cei,1:9,'air');
[b_ddata,b_tdata] = getDataContexts(cei,1:9,'belt');
[moo_ddata,moo_tdata] = getDataContexts(cei,1:9,'motionOnsetsOffsets');


data = {a_ddata{2} moo_ddata{2:4}};

sc = 0;

if sc > 0
    markerst = tei.contexts(sc).markers.air_onsets(1);
    markerse = tei.contexts(sc).markers.air_offsets(end);
    allFrames = b.frames_f(find(b.frames_f > markerst & b.frames_f < markerse));
else
    allFrames = b.frames_f;
end

n = 0;

speed = b.fSpeed;
if max(speed) > 1000
    speed(speed>35) = 35;
end
speed(speed<0) = 0;
% figure(101);clf;plot(speed)
minS = 0;
maxS = max(speed);

nBins = 30;
sBin = 3;
bins = linspace(minS,maxS,nBins);

ind = 1;
for ii = sBin:(length(bins)-1)
    if ii == 15
        n = 0;
    end
    st = bins(ii);
    se = bins(ii+1);
    thisB = (speed > st & speed < se);
    speedOnsets = find_rising_edge(thisB,0.5,1);
    speedOffsets = find_falling_edge(thisB,-0.5,1);
    [speedOnsets,speedOffsets] = resolveMotionOnsetsOffsets(speedOnsets,speedOffsets);
    [integrity pl] = checkIntegrity(speedOnsets,speedOffsets);
    if integrity < 0
        error
    end
    fns = [];
    dns = [];
    for jj = 1:length(speedOnsets)
        st1 = speedOnsets(jj);
        se1 = speedOffsets(jj);
        frames = find(allFrames > st1 & allFrames<se1);
        dur = b.ts(se1) - b.ts(st1);
        if length(frames) < 3
            continue;
        else
%             fns = [fns frames'];
            dns = [dns dur];
            fns{length(dns)} = frames;
        end
    end
    if isempty(fns)
        continue;
    end
    allFns{ind} = fns;
    allDns{ind} = dns;
    tbins(ind) = bins(ii);
    ind = ind + 1;
end
bins = tbins;
n=0;

signals = tei.deconv.spSigAll;
signals = [signals ei{2}.deconv.spSigAll];
inds = [];
for ii = 1:length(signals)
    thisSignal = signals{ii};
    for jj = 1:length(allFns)
        fns = allFns{jj};
        dns = allDns{jj};
        for kk = 1:length(fns)
            frames = fns{kk};
            spSig(kk) = mean(thisSignal(frames))/dns(kk);
%             semsptun(ii,jj) = std(thisSignal(frames))/sqrt(length(frames));
        end
        sptun(ii,jj) = mean(spSig);
    end
    thisSpeed = sptun(ii,:);
    [pks,locs,wpk,pkp] = findpeaks(thisSpeed,'MinPeakDistance',5,'MinPeakProminence',5);
    if length(pks) == 1
        inds = [inds ii];
%         tspeed(length(inds)) = bins(locs);
%         figure(101);clf;
%         plot(bins,thisSpeed);hold on;
%         plot(bins(locs),pks,'r*');
%         plotRastersMulti(data,ii,0,0,0);
%         n = 0;
    end
end
% inds = isnan(sptun);
% [rows,cols] = find(inds);
% sptun = sptun(:,(max(cols)+1):end);
% sBin = sBin + max(cols);
n = 0;

for ii = 1:length(signals)
%     figure(101);clf;
%     plot(bins(sBin:end),sptun(ii,:));hold on;
%     shadedErrorBar(bins(sBin:end),sptun(ii,:),semsptun(ii,:),{'color','b'},0.7);
%     errorbar(sptun(ii,:),semsptun(ii,:));
    MI(ii) = MutualInformation(bins',sptun(ii,:)');
%     [output,FRbin] = info_metrics_S(sptun(ii,:),bins(sBin:end),nBins-sBin+1,ones(size(sptun(ii,:))),100);
%     zMI(ii) = output.ShannonMI_Zsh;
end

n=0;

sptun = sptun(inds,:);

[~,peakPos] = max(sptun,[],2);
[~,cellNums] = sort(peakPos);
sptun = sptun(cellNums,:);
CRc = corrcoef(sptun);

figure(101);clf;
subplot(1,2,1)
imagesc(sptun,[0 max(sptun(:))/1.5]);colorbar;set(gca,'Ydir','Normal');
subplot(1,2,2)
imagesc(CRc);colorbar;set(gca,'Ydir','Normal');
% subplot(2,2,[3 4])
% plot(bins,sptun);





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function edge = find_rising_edge(signal,threshold,minimum_diff)
edge = find(diff(signal) >= threshold);
temp = find(diff(edge)<=minimum_diff);
edge(temp+1) = [];

function edge = find_falling_edge(signal,threshold,minimum_diff)
edge = find(diff(signal) <= threshold);
temp = find(diff(edge)<=minimum_diff);
edge(temp+1) = [];


% integriy check ... that is see if every motion onset has a corresponding motion offset
function [integrity ii] = checkIntegrity(speedOnsets,speedOffsets)
integrity = 1;
for ii = 1:(length(speedOnsets)-1)
    st = speedOnsets(ii); se = speedOnsets(ii+1);
    inds = find(speedOffsets > st & speedOffsets < se);
    if length(inds) > 1
        integrity = -1;
        return;
    end
end


for ii = 1:(length(speedOffsets)-1)
    st = speedOffsets(ii); se = speedOffsets(ii+1);
    inds = find(speedOnsets > st & speedOnsets < se);
    if length(inds) > 1
        integrity = -2;
        return;
    end
end


function [speedOnsets,speedOffsets] = resolveMotionOnsetsOffsets(speedOnsets,speedOffsets)
while 1 % see if there is a first motion offset which is less than the first motion onset
    if speedOffsets(1) < speedOnsets(1)
        speedOffsets(1) = [];
    else
        break;
    end
end
while 1 % see if there is an isolated last motion onset which is greater than the last motion onset
    if speedOnsets(end) > speedOffsets(end)
        speedOnsets(end) = [];
    else
        break;
    end
end