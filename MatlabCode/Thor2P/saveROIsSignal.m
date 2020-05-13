function saveROIsSignal (ei)

% clear all
% clc
% dataFolder = 'F:\Sam\Data\Animal_160513\CoverSlip940\05_29_2016\freeRun';
% ei = thorGetExperimentInfo(dataFolder);

if exist(ei.ROIsSignalFile)
    display('ROIs Signal File already exists');
    return;
end


ROIs = getROIPixels(ei);
rows = ei.pixelY;
cols = ei.pixelX;

hWaitBar = waitbar(0,sprintf('Processed Frames -'));
totalFrames = ei.totalFrames;
fid = fopen(ei.mcRawFile,'r');
ROIsSignal{length(ROIs)} = [];
for jj = 1:totalFrames
    waitbar(jj/totalFrames,hWaitBar,sprintf('Processing Frame %d',jj));
    fseek(fid, (jj-1)*rows*cols*2, 'bof');
    I1 = fread(fid,rows*cols,'uint16=>double',0,'l');
    I1 = reshape(I1,rows,cols)';
    for ii = 1:length(ROIs)
%         ROIsSignal{ii}.rectPixelValues(:,jj) = I1(ROIs{ii}.rectPixels);
        ROIsSignal.signalFromManualPixels(ii,jj) = mean(I1(ROIs{ii}.manualPixels));
    end
end

fclose(fid);
close(hWaitBar) ;

% windowSize = 
m = matfile(ei.ROIsSignalFile,'writable',true);
m.ROIsSignal = ROIsSignal;




% 
% nClusters = 2;
% ii = 1;
% while ii <= length(ROIs)
%     ii
%     thisROIPixels = ROIsSignal{ii}.pixelValues;
%     idx = kmeans(thisROIPixels,nClusters);
%     gr = unique(idx);
%     ROIsMask = zeros(size(averageImage));
% %     theOnes = ones(size(idx));
%     gri = 1;
%     while 1
% %         ROIsMask(ROIs{ii}(idx==gr(gri))) = theOnes(idx==gr(gri));
%         ROIsMask(ROIs{ii}(idx==gr(gri))) = 1;
%         maxI = max(averageImage(:));
%         figure(1);clf
%         subplot 121
%         sumImage = averageImage + ROIsMask*0.5*maxI;
%         axes('Position',[0.01 0.01 0.99 0.99]);
%         imagesc(sumImage);
%         colormap gray
%         axis off;
%         axis equal;
%         choice = questdlg('Verify','Verify Correctness of Pixels',...
%         'Yes','No','')
%         if strcmp(choice,'Yes')
%             ROIs{ii} = ROIs{ii}(idx == gr(gri));
%             nClusters = 2;
%             ii = ii + 1;
%             break;
%         else
%             gri = gri + 1;
%             if gri > length(gr);
%                 nClusters = nClusters + 1;
%                 break;
%             end
%         end
%     end
% end
% 
