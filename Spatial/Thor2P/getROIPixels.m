function ROIs = getROIPixels (expInfo)

% clear all
% clc
% 
% dataFolder = 'F:\Sam\Data\Animal_160513\CoverSlip940\05_29_2016\freeRun';
% expInfo = thorGetExperimentInfo(dataFolder);
[cvsROIs] = ReadImageJROI(expInfo.ROIsFile);

for ii = 1:length(cvsROIs)
    ROI = cvsROIs{ii};
    rect = ROI.vnRectBounds;
    boundary = ROI.mnCoordinates;
    [X,Y] = meshgrid(rect(2):rect(4),rect(1):rect(3));
    in = inpolygon(X,Y,boundary(:,1),boundary(:,2));
    rowIn = Y(in); colIn = X(in);
    X = reshape(X,numel(X),1);
    Y = reshape(Y,numel(Y),1);
    ROIs{ii}.rectPixels = sub2ind([expInfo.pixelY,expInfo.pixelX],Y,X);
    ROIs{ii}.manualPixels = sub2ind([expInfo.pixelY,expInfo.pixelX],rowIn,colIn);;
end

% ROIsMask = zeros(size(averageImage));
% theOnes = ones(size(colIn));
% indices = sub2ind(size(averageImage),rowIn,colIn);
% ROIsMask(indices) = theOnes;
% maxI = max(averageImage(:));
% figure(1);clf
% sumImage = averageImage + ROIsMask*0.5*maxI;
% imagesc(sumImage);
% colormap gray
% axis off;
% axis equal;