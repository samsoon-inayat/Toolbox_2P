function figureFrames(frames)

if ~exist('frames');
%     load('theFrames.mat');
    FileTif='Mean_G1_FL.tif';
    InfoImage=imfinfo(FileTif);
    mImage=InfoImage(1).Width;
    nImage=InfoImage(1).Height;
    NumberImages=length(InfoImage);

    FinalImage=zeros(nImage,mImage,NumberImages,'uint16');
    for i=1:NumberImages
       [FinalImage(:,:,i) map] =imread(FileTif,'Index',i);
    end
    frames = FinalImage;
%     frames = imread('Mean_G1_FL.tif',1:10);
end

figName = 'Fig'; thismfilename = mfilename; thismfilenamefullpath = mfilename('fullpath'); tmp = strfind(thismfilenamefullpath,'\');
figName = makeName(figName,thismfilenamefullpath(1:(tmp(end)-1)));
hf = figure(101);clf; pdfFile = 1; magFac = 1; columnWidth = magFac*6.8; columnHeight = magFac*1;
set(hf,'Position',[0.5 0.5 columnWidth columnHeight],'MenuBar','none','ToolBar','none',...
    'NumberTitle','on','Color','w','Units','inches','Resize','off',...
    'NextPlot','add');
fineSp = 0.001; coarseSp = 0.01; fontSize = 6; FillColorHists = 'k';

% frames = two_D_GaussianSourceTrajectorySink_two_fl();
% selectedFrames = linspace(1,numberOfFrames,numberOfFrames);%[1 2 3 300 305 310 350 355 356];
selectedFrames = ceil(linspace(30,60,17))
numberOfFrames = length(selectedFrames);

xs = linspace(0,0.85,numberOfFrames) - 0.03;
ys = ones(size(xs)) * 0.05;
height = 0.2;
width = 0.3;

% figure(1);clf;
for  ii = 1:length(selectedFrames)
    axes('Position',[xs(ii) ys(ii) height width]);
    imagesc(frames(:,:,selectedFrames(ii)));
    caxis([0 0.5]);
    box off;
    axis equal;
    axis off;
    colormap gray
end


fileName = sprintf('%s.pdf',figName);
if nameExists(fileName)
    recycle(fileName);
end
set(gcf, 'PaperPositionMode', 'Auto'); 
eval(sprintf('print -painters -dpdf -r600 ''%s''',fileName));
