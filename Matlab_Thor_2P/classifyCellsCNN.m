function pfs = classifyCellsCNN(allRs,ids,selCells)
plotFlag = 1;
Rs = allRs{ids(1)};
n = 0;

rows = size(Rs.rasters(:,:,1),1);
cols = size(Rs.rasters(:,:,1),2);

layers = [
    imageInputLayer([rows cols 1])
    
    convolution2dLayer(3,8,'Padding',1)
    batchNormalizationLayer
    reluLayer
    
    maxPooling2dLayer(2,'Stride',2)
    
    convolution2dLayer(3,16,'Padding',1)
    batchNormalizationLayer
    reluLayer
    
    maxPooling2dLayer(2,'Stride',2)
    
    convolution2dLayer(3,32,'Padding',1)
    batchNormalizationLayer
    reluLayer
    
    fullyConnectedLayer(10)
    softmaxLayer
    classificationLayer];

