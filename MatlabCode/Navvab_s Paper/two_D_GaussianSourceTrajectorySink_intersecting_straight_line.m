function two_D_GaussianSourceTrajectorySink_intersecting_straight_line
numberOfGridPoints = 128; % number of pixels in each x and y directions
frameRate = 150; % hz
frameCounter = 1;
durations = [0.2 0.2 0.1 0.1 0.25]; % durations of source, trajectory, and sink
sigmaSource = 5;

%% source pulse - duration = 0.5 sec
duration = durations(1);
numberOfFrames = round(duration * frameRate);
xo = -30; yo = 30;  % center of gaussian
xo1 = 30; yo1 = 30;  % center of gaussian
A = 1;
sigmaX = linspace(0.01,sigmaSource,numberOfFrames); sigmaY = sigmaX;
sigmaX1 = linspace(0.01,2*sigmaSource,numberOfFrames); sigmaY1 = sigmaX1;
gaussParam.center = [xo yo];
gaussParam.amplitude = A;
gaussParam.center1 = [xo1 yo1];
gaussParam.amplitude = A;
for ii = 1:length(sigmaX)
        gaussParam.sigma = [sigmaX(ii) sigmaY(ii)];
        gaussParam.sigma1 = [sigmaX1(ii) sigmaY1(ii)];
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
radius1 = sqrt(xo1^2 + yo1^2);
% startingAngle = -pi/2 + (-pi/4);
% startingAngle = atan2(yo,xo);
% endingAngle = startingAngle + pi/4;
% % startingAngle1 = pi/2 + (pi/4);
% startingAngle1 = atan2(yo1,xo1);
% endingAngle1 = startingAngle1 - pi/2;
% angles = linspace(startingAngle,endingAngle,numberOfFrames);
% angles1 = linspace(startingAngle1,endingAngle1,numberOfFrames);
gaussParam.amplitude = A;
gaussParam.sigma = [sigmaSource sigmaSource];
gaussParam.sigma1 = [2*sigmaSource 2*sigmaSource];
radii = linspace(radius,0,numberOfFrames);
radii1 = linspace(radius1,0,numberOfFrames);
for ii = 1:length(radii)
    thisAngle = pi/2+pi/4;%angles(ii);
    thisAngle1 = pi/4;%angles1(ii);
    xc = ceil(radii(ii) * cos(thisAngle));
    yc = ceil(radii(ii) * sin(thisAngle));
    xc1 = ceil(radii1(ii) * cos(thisAngle1));
    yc1 = ceil(radii1(ii) * sin(thisAngle1));
    gaussParam.center = [xc yc];
    gaussParam.center1 = [xc1 yc1];
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


%% Traveling Trajectory
duration = durations(2);
numberOfFrames = round(duration * frameRate);
radius = sqrt(xo^2 + yo^2);
radius1 = sqrt(xo1^2 + yo1^2);
% startingAngle = -pi/2 + (-pi/4);
% startingAngle = atan2(yo,xo);
% endingAngle = startingAngle + pi/4;
% % startingAngle1 = pi/2 + (pi/4);
% startingAngle1 = atan2(yo1,xo1);
% endingAngle1 = startingAngle1 - pi/2;
% angles = linspace(startingAngle,endingAngle,numberOfFrames);
% angles1 = linspace(startingAngle1,endingAngle1,numberOfFrames);
gaussParam.amplitude = A;
gaussParam.sigma = [sigmaSource sigmaSource];
gaussParam.sigma1 = [2*sigmaSource 2*sigmaSource];
radii = linspace(0,radius,numberOfFrames);
radii1 = linspace(0,radius1,numberOfFrames);
for ii = 1:length(radii)
    thisAngle = -pi/4;%angles(ii);
    thisAngle1 = pi+pi/4;%angles1(ii);
    xc = ceil(radii(ii) * cos(thisAngle));
    yc = ceil(radii(ii) * sin(thisAngle));
    xc1 = ceil(radii1(ii) * cos(thisAngle1));
    yc1 = ceil(radii1(ii) * sin(thisAngle1));
    gaussParam.center = [xc yc];
    gaussParam.center1 = [xc1 yc1];
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


%% Traveling Trajectory
duration = durations(2);
numberOfFrames = round(duration * frameRate);
radius = sqrt(xo^2 + yo^2);
radius1 = sqrt(xo1^2 + yo1^2);
% startingAngle = -pi/2 + (-pi/4);
startingAngle = atan2(yc,xc);
endingAngle = startingAngle + pi/4;
% startingAngle1 = pi/2 + (pi/4);
startingAngle1 = atan2(yc1,xc1);
endingAngle1 = startingAngle1 - pi/2;
angles = linspace(startingAngle,endingAngle,numberOfFrames);
angles1 = linspace(startingAngle1,endingAngle1,numberOfFrames);
gaussParam.amplitude = A;
gaussParam.sigma = [sigmaSource sigmaSource];
gaussParam.sigma1 = [2*sigmaSource 2*sigmaSource];
for ii = 1:length(angles)
    thisAngle = angles(ii);
    thisAngle1 = angles1(ii);
    xc = ceil(radius * cos(angles(1)));
    yc = ceil(radius * sin(thisAngle));
    xc1 = ceil(radius1 * cos(thisAngle1));
    yc1 = ceil(radius1 * sin(thisAngle1));
    gaussParam.center = [xc yc];
    gaussParam.center1 = [xc1 yc1];
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

%% Traveling Trajectory
duration = durations(3);
numberOfFrames = round(duration * frameRate);
xo = xc; yo = yc;
xo1 = xc1; yo1 = yc1;
radius = sqrt(xo^2 + yo^2);
radius1 = sqrt(xo1^2 + yo1^2);
startingAngle = endingAngle;
endingAngle = startingAngle + pi/4;
startingAngle1 = endingAngle1;
endingAngle1 = startingAngle1-pi/2;
angles = linspace(startingAngle,endingAngle,numberOfFrames);
angles1 = linspace(startingAngle1,endingAngle1,numberOfFrames);
gaussParam.amplitude = A;
gaussParam.sigma = [sigmaSource sigmaSource];
gaussParam.sigma1 = [2*sigmaSource 2*sigmaSource];
for ii = 1:length(angles)
    thisAngle = angles(ii);
    thisAngle1 = angles1(ii);
    xc = ceil(radius * cos(angles(1)));
    yc = ceil(radius * sin(thisAngle));
    xc1 = ceil(radius1 * cos(thisAngle1));
    yc1 = ceil(radius1 * sin(thisAngle1));
    gaussParam.center = [xc yc];
    gaussParam.center1 = [xc1 yc1];
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

%% Traveling Trajectory
duration = durations(4);
numberOfFrames = round(duration * frameRate);
radius = sqrt(xo^2 + yo^2);
radius1 = sqrt(xo1^2 + yo1^2);
startingAngle = endingAngle;
endingAngle = startingAngle + pi/2;
startingAngle1 = endingAngle1;
endingAngle1 = startingAngle1 - pi/2;
angles = linspace(startingAngle,endingAngle,numberOfFrames);
angles1 = linspace(startingAngle1,endingAngle1,numberOfFrames);
gaussParam.amplitude = A;
gaussParam.sigma = [sigmaSource sigmaSource];
gaussParam.sigma1 = [2*sigmaSource 2*sigmaSource];
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
    figure(1);clf;
    % surf(X,Y,gaussXY);
    imagesc(frames(:,:,frameCounter));
    colormap jet;
    axis off;
    axis equal;
    frameCounter = frameCounter + 1;
end

%% Sinking
duration = durations(5);
numberOfFrames = round(duration * frameRate);
xo = xc; yo = yc;  % center of gaussian
xo1 = xc1; yo1 = yc1;
A = 1;
sigmaX = linspace(sigmaSource,0.01,numberOfFrames); sigmaY = sigmaX;
sigmaX1 = linspace(2*sigmaSource,0.01,numberOfFrames); sigmaY1 = sigmaX1;
gaussParam.center = [xo yo];
gaussParam.center1 = [xo1 yo1];
gaussParam.amplitude = A;
for ii = 1:length(sigmaX)
        gaussParam.sigma = [sigmaX(ii) sigmaY(ii)];
        gaussParam.sigma1 = [sigmaX1(ii) sigmaY1(ii)];
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
xo1 = gaussParam.center1(1) ; yo1 = gaussParam.center1(2) ;
A = gaussParam.amplitude;
sigmaX = gaussParam.sigma(1); sigmaY = gaussParam.sigma(2);
sigmaX1 = gaussParam.sigma1(1); sigmaY1 = gaussParam.sigma1(2);
temp = linspace(ceil(-numberOfGridPoints/2)+1,ceil(numberOfGridPoints/2),numberOfGridPoints);
[X,Y] = meshgrid(temp,temp);
xTerm = ((X-xo).^2)./(2*(sigmaX^2));
yTerm = ((Y-yo).^2)./(2*(sigmaY^2));
xTerm1 = ((X-xo1).^2)./(2*(sigmaX1^2));
yTerm1 = ((Y-yo1).^2)./(2*(sigmaY1^2));
frame = A * exp(-(xTerm + yTerm)) + A * exp(-(xTerm1 + yTerm1));
