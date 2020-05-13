function two_D_GaussianSourceTrajectorySink_variable_Velocity
numberOfGridPoints = 128; % number of pixels in each x and y directions
frameRate = 150; % hz
frameCounter = 1;
durations = [0.25 1 0.5]; % durations of source, trajectory, and sink
sigmaSource = 5;

%% source pulse - duration = 0.5 sec
duration = durations(1);
numberOfFrames = round(duration * frameRate);
xo = -30; yo =-30;  % center of gaussian
A = 1;
sigmaX = linspace(0.01,sigmaSource,numberOfFrames); sigmaY = sigmaX;
gaussParam.center = [xo yo];
gaussParam.amplitude = A;
for ii = 1:length(sigmaX)
        gaussParam.sigma = [sigmaX(ii) sigmaY(ii)];
        gaussXY = generateFrame(numberOfGridPoints,gaussParam);
        frames(:,:,frameCounter) = flipud(gaussXY); 
        figure(1);clf;
        % surf(X,Y,gaussXY);
        imagesc(frames(:,:,frameCounter));
        colormap jet;
        axis off;
        axis equal;
        frameCounter = frameCounter + 1;
end

nothin = 0;

%% Traveling Trajectory
duration = durations(2);
numberOfFrames = round(duration * frameRate);
radius = sqrt(xo^2 + yo^2);
startingAngle = atan2(yo,xo);%-pi/2 + (-pi/4);
endingAngle = pi/2;
vf = 0.015;
tFrames = 0:numberOfFrames;
w = sin(2*pi*vf*tFrames);
figure(2);clf;
plot(tFrames,w)
theta = -cos(2*pi*vf*tFrames)+startingAngle + 1;
% angles = linspace(startingAngle,endingAngle,numberOfFrames);
angles = theta;
gaussParam.amplitude = A;
gaussParam.sigma = [sigmaSource sigmaSource];
for ii = 1:length(angles)
    thisAngle = angles(ii);
    xc = ceil(radius * cos(thisAngle));
    yc = ceil(radius * sin(thisAngle));
    gaussParam.center = [xc yc];
    gaussXY = generateFrame(numberOfGridPoints,gaussParam);
    frames(:,:,frameCounter) = flipud(gaussXY); 
    figure(1);clf;
    % surf(X,Y,gaussXY);
    imagesc(frames(:,:,frameCounter));
    colormap jet;
    axis off;
    axis equal;
    frameCounter = frameCounter + 1;
end

%% Sinking
duration = durations(3);
numberOfFrames = round(duration * frameRate);
xo = xc; yo = yc;  % center of gaussian
A = 1;
sigmaX = linspace(sigmaSource,0.01,numberOfFrames); sigmaY = sigmaX;
gaussParam.center = [xo yo];
gaussParam.amplitude = A;
for ii = 1:length(sigmaX)
        gaussParam.sigma = [sigmaX(ii) sigmaY(ii)];
        gaussXY = generateFrame(numberOfGridPoints,gaussParam);
        frames(:,:,frameCounter) = flipud(gaussXY); 
        figure(1);clf;
        % surf(X,Y,gaussXY);
        imagesc(frames(:,:,frameCounter));
        colormap jet;
        axis off;
        axis equal;
        frameCounter = frameCounter + 1;
end

frameCounter


function frame = generateFrame (numberOfGridPoints,gaussParam)
% this function generates a 2D gaussian provided the gauss param and number
% of grid points
xo = gaussParam.center(1) ; yo = gaussParam.center(2) ;
A = gaussParam.amplitude;
sigmaX = gaussParam.sigma(1); sigmaY = gaussParam.sigma(2);
temp = linspace(ceil(-numberOfGridPoints/2)+1,ceil(numberOfGridPoints/2),numberOfGridPoints);
[X,Y] = meshgrid(temp,temp);
xTerm = ((X-xo).^2)./(2*(sigmaX^2));
yTerm = ((Y-yo).^2)./(2*(sigmaY^2));
frame = A * exp(-(xTerm + yTerm));
