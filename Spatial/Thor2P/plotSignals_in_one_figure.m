function plotSignals_in_one_figure(ts,signal,varargin)

if nargin == 3
    fn = varargin{1};
else
    fn = 0;
end

oneX = ts;%1:size(signal,2);
numberOfROIs = size(signal,1);

xs = [];
ys = [];
xs = [xs oneX];
ys = [ys signal(1,:)/max(signal(1,:))];
xs = [xs NaN];
ys = [ys NaN];
figure(1+fn);clf;
for jj = 2:numberOfROIs
    xs = [xs oneX];
    oneSignal = signal(jj,:)/max(signal(jj,:));
    maxys = max(ys);
    ys = [ys (oneSignal + maxys + 1)];
    xs = [xs NaN];
    ys = [ys NaN];
end

figure(1+fn);clf;
plot(xs,ys,'LineWidth',0.1);

