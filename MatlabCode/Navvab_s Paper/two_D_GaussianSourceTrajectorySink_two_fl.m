function two_D_GaussianSourceTrajectorySink_two_fl ()

fn = 1000;

numberOfGridPoints = 128; % number of pixels in each x and y directions
frameRate = 150; % hz
frameCounter = 1;
durations = [0.25 0.25 1 0.5 0.5]; % durations of source, trajectory, and sink
sigmaSource = 5;

%% source pulse - duration = 0.5 sec
duration = durations(1);
numberOfFrames = round(duration * frameRate);
xo = -20; yo = 0;  % center of gaussian
xo1 = 20; yo1 = 0;  % center of gaussian
A = 1;
sigmaX = linspace(0.01,sigmaSource,numberOfFrames); sigmaY = sigmaX;
sigmaX1 = linspace(0.01,2*sigmaSource,numberOfFrames); sigmaY1 = sigmaX1;
gaussParam.center = [xo yo];
gaussParam.amplitude = A;
gaussParam.center1 = [xo1 yo1];
gaussParam.amplitude1 = A;
for ii = 1:length(sigmaX)
        gaussParam.sigma = [sigmaX(ii) sigmaY(ii)];
        gaussParam.sigma1 = [sigmaX1(ii) sigmaY1(ii)];
        gaussXY = generateFrame(numberOfGridPoints,gaussParam);
        frames(:,:,frameCounter) = flipud(gaussXY); 
        figure(fn);clf;
        % surf(X,Y,gaussXY);
        imagesc(frames(:,:,frameCounter));
        colormap jet;
        axis off;
        axis equal;
        title(sprintf('%d',frameCounter));
        frameCounter = frameCounter + 1;
end

nothin = 0;

%% Traveling Trajectory
duration = durations(2);
numberOfFrames = round(duration * frameRate);
radius = sqrt(xo^2 + yo^2);
radius1 = sqrt(xo1^2 + yo1^2);
startingAngle = atan2(yo,xo);
endingAngle = startingAngle + pi/3;
startingAngle1 = atan2(yo1,xo1);
endingAngle1 = startingAngle1 + pi/10;
angles = linspace(startingAngle,endingAngle,numberOfFrames);
angles1 = linspace(startingAngle1,endingAngle1,numberOfFrames);
gaussParam.amplitude = A;
sigmas = linspace(sigmaSource,3*sigmaSource,numberOfFrames); 
sigmas1 = linspace(2*sigmaSource,4*sigmaSource,numberOfFrames); 
for ii = 1:length(angles)
    thisAngle = angles(ii);
    thisAngle1 = angles1(ii);
    xc = ceil(radius * cos(thisAngle));
    yc = ceil(radius * sin(thisAngle));
    xc1 = ceil(radius1 * cos(thisAngle1));
    yc1 = ceil(radius1 * sin(thisAngle1));
    gaussParam.center = [xc yc];
    gaussParam.center1 = [xc1 yc1];
    gaussParam.sigma = [sigmaSource sigmas(ii)];
    gaussParam.sigma1 = [2*sigmaSource sigmas1(ii)];
    gaussXY = generateFrame(numberOfGridPoints,gaussParam);
    frames(:,:,frameCounter) = flipud(gaussXY); 
    figure(fn);clf;
    % surf(X,Y,gaussXY);
    imagesc(frames(:,:,frameCounter));
    colormap jet;
    axis off;
    axis equal;
    title(sprintf('%d',frameCounter));
    frameCounter = frameCounter + 1;
end

%% Traveling Trajectory
duration = durations(3);
numberOfFrames = round(duration * frameRate);
radius = sqrt(xo^2 + yo^2);
radius1 = sqrt(xo1^2 + yo1^2);
startingAngle = endingAngle;
endingAngle = startingAngle + pi/3;
startingAngle1 = endingAngle1;
endingAngle1 = startingAngle1+pi/10;
angles = linspace(startingAngle,endingAngle,numberOfFrames);
angles1 = linspace(startingAngle1,endingAngle1,numberOfFrames);
gaussParam.amplitude = A;
gaussParam.sigma = [sigmaSource sigmas(end)];
gaussParam.sigma1 = [2*sigmaSource sigmas1(end)];
for ii = 1:length(angles)
    thisAngle = angles(ii);
    thisAngle1 = angles1(ii);
    xc = ceil(radius * cos(thisAngle));
    yc = ceil(radius * sin(thisAngle));
    xc1 = ceil(radius1 * cos(thisAngle1));
    yc1 = ceil(radius1 * sin(thisAngle1));
    gaussParam.center = [xc yc];
    gaussParam.center1 = [xc1 yc1];
    gaussXY = generateFrame(numberOfGridPoints,gaussParam);
    frames(:,:,frameCounter) = flipud(gaussXY); 
    figure(fn);clf;
    % surf(X,Y,gaussXY);
    imagesc(frames(:,:,frameCounter));
    colormap jet;
    axis off;
    axis equal;
    title(sprintf('%d',frameCounter));
    frameCounter = frameCounter + 1;
end

%% Traveling Trajectory
duration = durations(4);
numberOfFrames = round(duration * frameRate);
radius = sqrt(xo^2 + yo^2);
radius1 = sqrt(xo1^2 + yo1^2);
startingAngle = endingAngle;
endingAngle = startingAngle + pi/2;
startingAngle1 = endingAngle1;
endingAngle1 = startingAngle1 + pi/10;
angles = linspace(startingAngle,endingAngle,numberOfFrames);
angles1 = linspace(startingAngle1,endingAngle1,numberOfFrames);
gaussParam.amplitude = A;
% gaussParam.sigma = [sigmaSource sigmaSource];
% gaussParam.sigma1 = [2*sigmaSource 2*sigmaSource];
for ii = 1:length(angles)
    thisAngle = angles(ii);
    thisAngle1 = angles1(ii);
    xc = ceil(radius * cos(thisAngle));
    yc = ceil(radius * sin(thisAngle));
    xc1 = ceil(radius1 * cos(thisAngle1));
    yc1 = ceil(radius1 * sin(thisAngle1));
    gaussParam.center = [xc yc];
    gaussParam.center1 = [xc1 yc1];
    gaussXY = generateFrame(numberOfGridPoints,gaussParam);
    frames(:,:,frameCounter) = flipud(gaussXY); 
    figure(fn);clf;
    % surf(X,Y,gaussXY);
    imagesc(frames(:,:,frameCounter));
    colormap jet;
    axis off;
    axis equal;
    title(sprintf('%d',frameCounter));
    frameCounter = frameCounter + 1;
end

%% Sinking
duration = durations(5);
numberOfFrames = round(duration * frameRate);
xo = xc; yo = yc;  % center of gaussian
xo1 = xc1; yo1 = yc1;
A = 1;
sigmaX = linspace(sigmas(end),0.01,numberOfFrames); sigmaY = linspace(sigmaSource,0.01,numberOfFrames);;
sigmaX1 = linspace(2*sigmaSource,0.01,numberOfFrames); sigmaY1 = sigmaX1;
gaussParam.center = [xo yo];
gaussParam.center1 = [xo1 yo1];
gaussParam.amplitude = A;
for ii = 1:length(sigmaX)
        gaussParam.sigma = [sigmaY(ii) sigmaX(ii)];
        gaussParam.sigma1 = [sigmaX1(ii) sigmaY1(ii)];
        gaussXY = generateFrame(numberOfGridPoints,gaussParam);
        frames(:,:,frameCounter) = flipud(gaussXY); 
        figure(fn);clf;
        % surf(X,Y,gaussXY);
        imagesc(frames(:,:,frameCounter));
        colormap jet;
        axis off;
        axis equal;
        title(sprintf('%d',frameCounter));
        frameCounter = frameCounter + 1;
end

frameCounter-1;


function frame = generateFrame (numberOfGridPoints,gaussParam)
% this function generates a 2D gaussian provided the gauss param and number
% of grid points
xo = gaussParam.center(1) ; yo = gaussParam.center(2) ;
xo1 = gaussParam.center1(1) ; yo1 = gaussParam.center1(2) ;
A = gaussParam.amplitude;
A1 = gaussParam.amplitude1;
sigmaX = gaussParam.sigma(1); sigmaY = gaussParam.sigma(2);
sigmaX1 = gaussParam.sigma1(1); sigmaY1 = gaussParam.sigma1(2);
temp = linspace(ceil(-numberOfGridPoints/2)+1,ceil(numberOfGridPoints/2),numberOfGridPoints);
[X,Y] = meshgrid(temp,temp);
xTerm = ((X-xo).^2)./(2*(sigmaX^2));
yTerm = ((Y-yo).^2)./(2*(sigmaY^2));
xTerm1 = ((X-xo1).^2)./(2*(sigmaX1^2));
yTerm1 = ((Y-yo1).^2)./(2*(sigmaY1^2));
frame = A * exp(-(xTerm + yTerm)) + A1* exp(-(xTerm1 + yTerm1));
frame = frame/max(frame(:));
