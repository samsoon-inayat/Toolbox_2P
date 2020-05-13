function plotROIs_all_in_one_figure(ei,varargin)

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
    ROIsSignal.signalFromManualPixels = signalFromManualPixels;
    m = matfile(ei.ROIsSignalFile,'Writable',true);
    m.ROIsSignal = ROIsSignal;
end
% figure(1);clf;
ROIs = getROIPixels(ei);
load(ei.averageImageFile);
numberOfROIs = length(ROIs);
meanSignal = mean(Image_0001_0001_MC_averageImage(:));
signal = (signalFromManualPixels - meanSignal)/meanSignal;

oneX = 1:size(signal,2);

xs = [];
ys = [];
xs = [xs oneX];
ys = [ys signal(1,:)/max(signal(1,:))];
figure(1+fn);clf;
for jj = 2:numberOfROIs
    xs = [xs oneX];
    oneSignal = signal(jj,:)/max(signal(jj,:));
    maxys = max(ys);
    ys = [ys (oneSignal + maxys + 1)];
end

figure(1+fn);clf;
plot(xs,ys,'LineWidth',0.1);

