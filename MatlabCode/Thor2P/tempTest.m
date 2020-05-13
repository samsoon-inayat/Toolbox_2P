
% m = matfile(ei.ROIsSignalFile);
% ROIsSignal = m.ROIsSignal;
% figure(1);clf;
% ROIs = getROIPixels(ei);
% load(ei.averageImageFile);
% numberOfROIs = length(ROIs);
% left = 0.03;
% width = 0.9;
% height = 0.9/25;
% bottoms = linspace(1-height,0,25);
% verticalSpace = 0.001;
% meanSignal = mean(Image_0001_0001_MC_averageImage(:));
% for ii = 1:25
%     axes('Position',[left bottoms(ii) width height]);
%     plot((ROIsSignal{ii}.signalFromManualPixels-meanSignal)/meanSignal); hold on;
% %     meanRectSignal = mean(ROIsSignal{ii}.rectPixelValues);
% %     plot((meanRectSignal - meanSignal)/meanSignal,'r');
%     set(gca,'xtick',[],'ytick',[],'box','off','ylim',[-0.5 11]);
%     ylabel(sprintf('%d',ii),'FontSize',8);
%     
% %     if ii < length(ROIs)
% %         set(gca,'xtick',[]);
% %     end
% end
% % 
% % figure(1);clf;
% % ii = 20;
% % plot((ROIsSignal{ii}.signalFromManualPixels-meanSignal)/meanSignal); hold on;
% % meanRectSignal = mean(ROIsSignal{ii}.rectPixelValues);
% % plot((meanRectSignal - meanSignal)/meanSignal,'r');
% 
% % 
% % Image_0001_0001_MC_maxImage_1 = reshape(Image_0001_0001_MC_maxImage,ei.pixelY,ei.pixelX);
% % Image_0001_0001_MC_maxImage_1 = Image_0001_0001_MC_maxImage_1';
% % figure(1);clf;
% % imagesc(Image_0001_0001_MC_maxImage_1);
% % axis off;
% % axis equal;
% % m1 = matfile(ei.averageImageFile);
% % m1.Image_0001_0001_MC_maxImage = Image_0001_0001_MC_maxImage_1;
% % imwrite(uint16(Image_0001_0001_MC_maxImage_1),ei.maxImageTifFile,'tif');