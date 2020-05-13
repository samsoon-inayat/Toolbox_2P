function two_D_GaussianSourceTrajectorySink_multiple

numberOfGridPoints = 128; % number of pixels in each x and y directions
frameRate = 150; % hz
frameCounter = 1;
durations{1} = [0.25 1 0.1];
durations{2} =  [0.1 1 0.25];
startTimeOfSources = [0 0.35];
durations1 = sum(durations{1});
durations2 = sum(durations{2});
totalTime = max(durations1 + startTimeOfSources(1),durations2 + startTimeOfSources(2));
% totalTime = durations1;
totalFrames =round( totalTime * frameRate);
startFrameOfSources = round(startTimeOfSources * frameRate) + 1;
numberOfSources = 2;

for ii = 1:numberOfSources
    numberOfFrames{ii} = round(durations{ii} * frameRate);
end
for ii = 1:numberOfSources
    totalFramesSources(ii) = sum(numberOfFrames{ii});
end

visualization = [0 0];
gP = source1(numberOfFrames{1});
frames{1} = generateFrames(numberOfGridPoints,gP,visualization(1));
gP = source2(numberOfFrames{2});
frames{2} = generateFrames(numberOfGridPoints,gP,visualization(2));


allFrames = frames{1};
zP = zeros([size(allFrames(:,:,1)) totalFrames-totalFramesSources(1)]);
allFrames1 = cat(3,allFrames,zP);

allFrames = frames{2};
zP = zeros([size(allFrames(:,:,1)) totalFrames-totalFramesSources(2)]);
allFrames2 = cat(3,zP,allFrames);

allFrames = allFrames1+allFrames2;
for ii = 1:totalFrames
    figure(1);clf;
    imagesc(allFrames(:,:,ii));
    axis off;
    axis equal;
    title(sprintf('%d',ii));
end

function frameSeq = generateFrames (numberOfGridPoints,gP,v)
allxo = gP(1,:);
allyo = gP(2,:);
allA = gP(3,:);
allsigmaX = gP(4,:);
allsigmaY = gP(5,:);
allRA = gP(6,:);

for ii = 1:length(allxo)
    gaussParam.center = [allxo(ii) allyo(ii)];
    gaussParam.amplitude = allA(ii);
    gaussParam.sigma = [allsigmaX(ii) allsigmaY(ii)];
    frameSeq(:,:,ii) = generateRotatedFrame(numberOfGridPoints,gaussParam,allRA(ii));
    if v
        figure(1);clf;
        imagesc(frameSeq(:,:,ii));
        axis equal
        axis off
    end
end

function frame = generateRotatedFrame (numberOfGridPoints,gaussParam,gRA)
% this function generates a 2D gaussian provided the gauss param and number
% of grid points
x0 = gaussParam.center(1) ; y0 = gaussParam.center(2) ;
A = gaussParam.amplitude;
sigma_x = gaussParam.sigma(1); sigma_y = gaussParam.sigma(2);
temp = linspace(ceil(-numberOfGridPoints/2)+1,ceil(numberOfGridPoints/2),numberOfGridPoints);
[X,Y] = meshgrid(temp,temp);
a = cos(gRA)^2/(2*sigma_x^2) + sin(gRA)^2/(2*sigma_y^2);
b = -sin(2*gRA)/(4*sigma_x^2) + sin(2*gRA)/(4*sigma_y^2);
c = sin(gRA)^2/(2*sigma_x^2) + cos(gRA)^2/(2*sigma_y^2);

frame = A*exp( - (a*(X-x0).^2 - 2*b*(X-x0).*(Y-y0) + c*(Y-y0).^2)) ;
frame  = flipud(frame);


function val = source1(numberOfFrames)
sigmaSource = [5];
sM =[2.5] ;
allxo = -20 * ones(1,numberOfFrames(1)); 
allyo =-20*ones(1,numberOfFrames(1));
allA = 1*ones(1,numberOfFrames(1));
allsigmaX = linspace(0.01,sigmaSource(1),numberOfFrames(1)); 
allsigmaY = sM(1)*allsigmaX;
allRA = zeros(1,numberOfFrames(1));
xo = allxo(end);
yo = allyo(end);
radius = sqrt(xo^2 + yo^2);
startingAngle = atan2(yo,xo);%-pi/2 + (-pi/4);
vf = 0.5/numberOfFrames(2);
tFrames = linspace(0,numberOfFrames(2),numberOfFrames(2));
w = sin(2*pi*vf*tFrames);
theta = -cos(2*pi*vf*tFrames)+startingAngle + 1;
angles = theta;
A = 1;
RA = linspace(0,2*pi,numberOfFrames(2));
for ii = 1:length(angles)
    thisAngle = angles(ii);
    xc = ceil(radius * cos(thisAngle));
    yc = ceil(radius * sin(thisAngle));
    allxo = [allxo xc];
    allyo = [allyo yc];
    allA = [allA A];
    allsigmaX = [allsigmaX sigmaSource(1)];
    allsigmaY = [allsigmaY sM(1)*sigmaSource(1)];
    allRA = [allRA RA(ii)];
end
xo = allxo(end); yo = allyo(end);  % center of gaussian
A = 1;
sigmaX = linspace(sigmaSource(1),0.01,numberOfFrames(3));
sigmaY = sM(1)*sigmaX;
for ii = 1:length(sigmaX)
    allxo = [allxo xc];
    allyo = [allyo yc];
    allA = [allA A];
    allsigmaX = [allsigmaX sigmaX(ii)];
    allsigmaY = [allsigmaY sigmaY(ii)];
    allRA = [allRA RA(end)];
end
val = [allxo;allyo;allA;allsigmaX;allsigmaY;allRA];

function val = source2(numberOfFrames)
sigmaSource = 6;
sM = 5;
allxo = 10 * ones(1,numberOfFrames(1)); 
allyo =10*ones(1,numberOfFrames(1));
allA = 1*ones(1,numberOfFrames(1));
allsigmaX = linspace(0.01,sigmaSource(1),numberOfFrames(1)); 
allsigmaY = sM(1)*allsigmaX;
allRA = zeros(1,numberOfFrames(1));
xo = allxo(end);
yo = allyo(end);
radius = sqrt(xo^2 + yo^2);
startingAngle = atan2(yo,xo);%-pi/2 + (-pi/4);
vf = 0.5/numberOfFrames(2);
tFrames = linspace(0,numberOfFrames(2),numberOfFrames(2));
w = sin(2*pi*vf*tFrames);
theta = -cos(2*pi*vf*tFrames)+startingAngle + 1;
angles = theta;
A = 1;
RA = linspace(0,2*pi,numberOfFrames(2));
for ii = 1:length(angles)
    thisAngle = angles(ii);
    xc = ceil(radius * cos(thisAngle));
    yc = ceil(radius * sin(thisAngle));
    allxo = [allxo xc];
    allyo = [allyo yc];
    allA = [allA A];
    allsigmaX = [allsigmaX sigmaSource(1)];
    allsigmaY = [allsigmaY sM(1)*sigmaSource(1)];
    allRA = [allRA RA(ii)];
end
xo = allxo(end); yo = allyo(end);  % center of gaussian
A = 1;
sigmaX = linspace(sigmaSource(1),0.01,numberOfFrames(3));
sigmaY = sM(1)*sigmaX;
for ii = 1:length(sigmaX)
    allxo = [allxo xc];
    allyo = [allyo yc];
    allA = [allA A];
    allsigmaX = [allsigmaX sigmaX(ii)];
    allsigmaY = [allsigmaY sigmaY(ii)];
    allRA = [allRA RA(end)];
end
val = [allxo;allyo;allA;allsigmaX;allsigmaY;allRA];