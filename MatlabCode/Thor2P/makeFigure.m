function makeFigure (abfd)

frameTrigger = abfd(:,1);
encoderA = abfd(:,2);
encoderB = abfd(:,3);
encoderR = abfd(:,4);
reward = abfd(:,5);
timeAxis = abfd(:,6);

figName = 'Fig'; thismfilename = mfilename; thismfilenamefullpath = mfilename('fullpath'); tmp = strfind(thismfilenamefullpath,'\');
figName = makeName(figName,thismfilenamefullpath(1:(tmp(end)-1)));
hf = figure(101);clf; pdfFile = 1; magFac = 1; columnWidth = magFac*6.8; columnHeight = magFac*10;
set(hf,'Position',[0.5 0.5 columnWidth columnHeight],'MenuBar','none','ToolBar','none',...
    'NumberTitle','on','Color','w','Units','inches','Resize','off',...
    'NextPlot','add');
fineSp = 0.001; coarseSp = 0.01; fontSize = 6; FillColorHists = 'k';

panelLabels = {'A' 'B' 'C' 'D'};
panelLefts = [0 0 0 0]+6*coarseSp; panelBottoms = [0.7 0.45 0.2 0]+5*coarseSp;
panelWidths = [1 1 1 1]*0.99; panelHeights = [1 1 1 1]*0.15;

% generating axes
pn = 1;
pLeft = panelLefts(pn); pBottom = panelBottoms(pn); pWidth = panelWidths(pn); pHeight = panelHeights(pn);
annotation('textbox',[pLeft-5*coarseSp,pBottom+pHeight+2*coarseSp,0.1,0.1],'String',panelLabels{pn},'EdgeColor','none','FontWeight','Bold');
sPanelLefts = [0];sPanelBottoms = [0];sPanelWidths = [0.8];sPanelHeights = [1];
sPanelLefts = pLeft + sPanelLefts * pWidth;sPanelBottoms = pBottom + sPanelBottoms * pHeight;
sPanelHeights = sPanelHeights * pHeight;sPanelWidths = sPanelWidths * pWidth;
ii = 1; left = sPanelLefts(ii); bottom = sPanelBottoms(ii) ; width = sPanelWidths(ii); height = sPanelHeights(ii);
axes('Position',[left bottom width height]);

plot(timeAxis,encoderR);

% generating axes
pn = 2;
pLeft = panelLefts(pn); pBottom = panelBottoms(pn); pWidth = panelWidths(pn); pHeight = panelHeights(pn);
annotation('textbox',[pLeft-5*coarseSp,pBottom+pHeight+2*coarseSp,0.1,0.1],'String',panelLabels{pn},'EdgeColor','none','FontWeight','Bold');
sPanelLefts = [0];sPanelBottoms = [0];sPanelWidths = [0.8];sPanelHeights = [1];
sPanelLefts = pLeft + sPanelLefts * pWidth;sPanelBottoms = pBottom + sPanelBottoms * pHeight;
sPanelHeights = sPanelHeights * pHeight;sPanelWidths = sPanelWidths * pWidth;
ii = 1; left = sPanelLefts(ii); bottom = sPanelBottoms(ii) ; width = sPanelWidths(ii); height = sPanelHeights(ii);
axes('Position',[left bottom width height]);


plot(timeAxis,encoderA);

pn = 3;
pLeft = panelLefts(pn); pBottom = panelBottoms(pn); pWidth = panelWidths(pn); pHeight = panelHeights(pn);
annotation('textbox',[pLeft-5*coarseSp,pBottom+pHeight+2*coarseSp,0.1,0.1],'String',panelLabels{pn},'EdgeColor','none','FontWeight','Bold');
sPanelLefts = [0];sPanelBottoms = [0];sPanelWidths = [0.8];sPanelHeights = [1];
sPanelLefts = pLeft + sPanelLefts * pWidth;sPanelBottoms = pBottom + sPanelBottoms * pHeight;
sPanelHeights = sPanelHeights * pHeight;sPanelWidths = sPanelWidths * pWidth;
ii = 1; left = sPanelLefts(ii); bottom = sPanelBottoms(ii) ; width = sPanelWidths(ii); height = sPanelHeights(ii);
axes('Position',[left bottom width height]);

plot(timeAxis,encoderB);


% generating axes
pn = 4;
pLeft = panelLefts(pn); pBottom = panelBottoms(pn); pWidth = panelWidths(pn); pHeight = panelHeights(pn);
annotation('textbox',[pLeft-5*coarseSp,pBottom+pHeight+2*coarseSp,0.1,0.1],'String',panelLabels{pn},'EdgeColor','none','FontWeight','Bold');
sPanelLefts = [0];sPanelBottoms = [0];sPanelWidths = [0.8];sPanelHeights = [1];
sPanelLefts = pLeft + sPanelLefts * pWidth;sPanelBottoms = pBottom + sPanelBottoms * pHeight;
sPanelHeights = sPanelHeights * pHeight;sPanelWidths = sPanelWidths * pWidth;
ii = 1; left = sPanelLefts(ii); bottom = sPanelBottoms(ii) ; width = sPanelWidths(ii); height = sPanelHeights(ii);
axes('Position',[left bottom width height]);

plot(timeAxis,reward);

% Printing
if pdfFile
    fileName = sprintf('%s.pdf',figName);
    if nameExists(fileName)
        recycle(fileName);
    end
    set(gcf, 'PaperPositionMode', 'Auto'); eval(sprintf('print -painters -dpdf -r600 ''%s''',fileName));
end