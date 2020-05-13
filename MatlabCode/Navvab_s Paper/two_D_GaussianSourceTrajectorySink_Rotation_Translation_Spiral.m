function two_D_GaussianSourceTrajectorySink_Rotation_Translation_Spiral

numberOfGridPoints = 128; % number of pixels in each x and y directions
frameRate = 150; % hz
frameCounter = 1;
durations = [0.1 0.5 0.1]; % durations of source, trajectory, and sink
sigmaSource = 5;
sM =2.5 ;

%% source pulse - duration = 0.5 sec
duration = durations(1);
numberOfFrames = round(duration * frameRate);
xo = -20; yo =-20;  % center of gaussian
A = 1;
sigmaX = linspace(0.01,sigmaSource,numberOfFrames); sigmaY = sM*sigmaX;
gaussParam.center = [xo yo];
gaussParam.amplitude = A;
for ii = 1:length(sigmaX)
        gaussParam.sigma = [sigmaX(ii) sigmaY(ii)];
        gaussXY = generateRotatedFrame(numberOfGridPoints,gaussParam,0);
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
vf = 0.5/numberOfFrames;
tFrames = linspace(0,numberOfFrames,numberOfFrames);
w = sin(2*pi*vf*tFrames);
figure(2);clf;
plot(tFrames,w)
theta = -cos(2*pi*vf*tFrames)+startingAngle + 1;
% angles = linspace(startingAngle,endingAngle,numberOfFrames);
angles = theta;
gaussParam.amplitude = A;
gaussParam.sigma = [sigmaSource sM*sigmaSource];
RA = linspace(0,2*pi,numberOfFrames);
for ii = 1:length(angles)
    thisAngle = angles(ii);
    xc = ceil(radius * cos(thisAngle));
    yc = ceil(radius * sin(thisAngle));
    gaussParam.center = [xc yc];
    gaussXY = generateRotatedFrame(numberOfGridPoints,gaussParam,RA(ii));
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
sigmaX = linspace(sigmaSource,0.01,numberOfFrames); sigmaY = sM*sigmaX;
gaussParam.center = [xo yo];
gaussParam.amplitude = A;
for ii = 1:length(sigmaX)
        gaussParam.sigma = [sigmaX(ii) sigmaY(ii)];
        gaussXY = generateRotatedFrame(numberOfGridPoints,gaussParam,RA(end));
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
%     figure(1);clf;
%     imagesc(Z)