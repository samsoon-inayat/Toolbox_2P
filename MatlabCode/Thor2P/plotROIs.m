function plotROIs(ei,varargin)

if nargin == 2
    fn = varargin{1};
else
    fn = 0;
end

m = matfile(ei.ROIsSignalFile);
ROIsSignal = m.ROIsSignal;
try
signalFromManualPixels = ROIsSignal.signalFromManualPixels;
catch
    for ii = 1:length(ROIsSignal)
        signalFromManualPixels(ii,:) = ROIsSignal{ii}.signalFromManualPixels;
    end
end
% figure(1);clf;
ROIs = getROIPixels(ei);
load(ei.averageImageFile);
numberOfROIs = length(ROIs);

totalNumber = numberOfROIs;
perSet = 20;
% numberOfSets = floor(totalNumber/perSet); % find how many sets of frames to be processed  with par for
startOfSet = 1:perSet:totalNumber;
endOfSet = (1:perSet:totalNumber)-1;
endOfSet(1) = []; endOfSet(end+1) = totalNumber;
left = 0.05;
width = 0.9;
height = 0.9/perSet;
bottoms = linspace(1-height,0,perSet);
verticalSpace = 0.001;
meanSignal = mean(Image_0001_0001_MC_averageImage(:));

numberOfSets = length(startOfSet);

for jj = 1:numberOfSets
    figure(jj+fn);clf;
    set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.
    roins = startOfSet(jj):endOfSet(jj);
    for iii = 1:perSet
        axes('Position',[left bottoms(iii) width height]);
        ii = roins(iii);
        signal = (signalFromManualPixels(ii,:)-meanSignal)/meanSignal;
        plot(signal); hold on;
        maxS = max(signal);
        minS = min(signal);
        set(gca,'xtick',[],'ytick',[],'box','off','ylim',[minS maxS]);
        ylabel(sprintf('%d',ii),'FontSize',8);
        text(floor(length(signal)/2),floor((maxS-minS)/2),sprintf('%.1f %.1f',minS,maxS));
    end
end
