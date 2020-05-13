function baselines = findBleachingTrend (signal,varargin)
nFrames = size(signal,1);
baselines = zeros(size(signal));
if ~isempty(varargin)
    pDegree = varargin{1};
else
    pDegree = 2;
end

% d1 = designfilt('lowpassiir','FilterOrder',12, ...
%     'HalfPowerFrequency',0.15,'DesignMethod','butter');

parfor ii = 1:size(signal,2)
    sig = signal(:,ii);
    if length(unique(sig)) == 1
        baselines(:,ii) = ones(nFrames,1) * sig(1);
        continue;
    end
    minimas1 = detectPeaksN(sig);
    if isempty(minimas1)
        baselines(:,ii) = ones(nFrames,1) * min(sig);
        continue;
    end
    x = minimas1(:,1);
    if length(x) == 1
        baselines(:,ii) = ones(nFrames,1) * sig(x);
        continue;
    end
%     y = minimas1(:,2);
    sigI = interp1(x,sig(x),min(x):max(x));
    ssigI = ones(size(1:(min(x)-1)))*sigI(1);
    esigI = ones(size((max(x)+1):nFrames))*sigI(end);
    finalSigI = [ssigI sigI esigI];
%     finalSigI = filtfilt(d1,finalSigI);
    p = polyfit(1:nFrames,finalSigI,pDegree);
%     p = polyfit(x,y,2);
    finalSigI = polyval(p,1:nFrames);
    baselines(:,ii) = finalSigI';
end


