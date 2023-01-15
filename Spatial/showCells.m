function ha = showCells(ha,tei,pl,selCells,perc)
%%
% if isempty(ha)
%     tei = evalin('base','tei(4)');
%     pl = 1;
% end
% tei = ei{4};
[ops,astat,arecells] = get_ops(tei,pl);
perc = [0.3 0.3];
selCells = logical(arecells(:,1));
% astat = astat(logical(arecells(:,1)));
FS = 8;
micronsPerPixel = tei.thorExp.widthUM/tei.thorExp.pixelX;
xrange = ops.xrange;
yrange = ops.yrange;
mimg = ops.meanImgE;
% mimg = ops.refImg;
if isempty(ha)
    figure(102);clf;plot(0,0);
    ha = gca;
end
figure(1000);clf;
    ha = axes;
axes(ha);
imagesc(mimg,[min(mimg(:)) max(mimg(:))]);
hold on;
maskZ = zeros(size(mimg));
colormap gray;
axis equal;
axis off;
hold on;
allmask = maskZ;
iii = 0;
for ii = 58%:length(astat)
  if ~logical(arecells(ii,1))
    continue;
  end
  iii = iii + 1;
  if ~selCells(iii)
    continue;
  end
  stat = astat(ii);
  mask = maskZ;
  mask(stat{1}.ipix) = 1;
  mask = mask';
  mask = expandOrCompressMask(mask,stat{1}.footprint);
  allmask(mask==1) = 1;
%   cellX = double(stat{1}.xpix);% + double(min(xrange));
%   cellY = double(stat{1}.ypix);% + double(min(yrange));
%   plot(cellX,cellY);
end
img3 = cat(3,mimg,mimg,mimg);
ch = 2;
rch = img3(:,:,ch);
rch = rch + max(rch(:))*(perc(1)*allmask);
img3(:,:,ch) = rch;
img2 = img3(:,:,setdiff([1 2 3],ch));
imagesc(img3,[min(img2(:)) 	perc(2)*max(img2(:))]);
% mx = min(cellX);
% my = max(cellY);
% text(mx-10,my+20,num2str(217),'color','c','fontsize',FS,'FontWeight','Bold');  
return;



if ~iscell(selCells)
    ccs = (selCells);%tei.plane{pl}.tP.areCells(logical(selCells));
    for cc = 1:length(ccs)
        if ~ccs(cc)
          continue;
        end
        stat =  astat(ccs(cc));
        mask = maskZ;
        mask(stat{1}.ipix) = 1;
%         mask = expandOrCompressMask(mask,stat{1}.footprint);
        mask = mask';
        allmask(mask==1) = 1;
    end
    imagesc(allmask);
%     img3 = cat(3,mimg,mimg,mimg);
%     ch = 2;
%     rch = img3(:,:,ch);
%     rch = rch + max(rch(:))*(perc(1)*allmask);
%     img3(:,:,ch) = rch;
%     img2 = img3(:,:,setdiff([1 2 3],ch));
%     imagesc(img3,[min(img2(:)) 	perc(2)*max(img2(:))]);hold on;
    for cc = 1:length(ccs)
        if ~ccs(cc)
          continue;
        end
        stat =  astat(ccs(cc));
        cellX = double(stat{1}.xpix);% + double(min(xrange));
        cellY = double(stat{1}.ypix);% + double(min(yrange));
        mx = min(cellX);
        my = max(cellY);
%         if selCells(cc) == 46
%             text(mx+30,my+0,num2str(selCells(cc)),'color','c','fontsize',FS+2,'FontWeight','Normal');
%         else
%             if selCells(cc) == 22
%                 text(mx-60,my-20,num2str(selCells(cc)),'color','c','fontsize',FS+2,'FontWeight','Normal');  
%             else
                text(mx-10,my+20,num2str(selCells(cc)),'color','c','fontsize',FS,'FontWeight','Bold');  
%             end
%         end
%         scatter(cellX,cellY,'.r');
%         plot(cellX,cellY,'.','color','r');
    end
else
    imagesc(mimg);hold on;
    accs = selCells; %clear selCells;
    cColors = {'c','r','g','b','m','y','w'};
%     cColors = {'c','r','g','r','g','y'};
    tcColors = {'c','r','g','b','m','y','w'};
    for ss = 1:length(accs)
        ccs = accs{ss};
        for cc = 1:length(ccs)
            stat =  astat(ccs(cc));
            cellX = double(stat{1}.xpix);% + double(min(xrange));
            cellY = double(stat{1}.ypix);% + double(min(yrange));
%             cellX = astat(ccs(cc)).xpix + min(xrange);
%             cellY = astat(ccs(cc)).ypix + min(yrange);
%             mx = min(cellX);
%             my = max(cellY);
%             text(mx-5,my,num2str(selCells(cc)),'color',tcColors{ss},'fontsize',FS);
            plot(cellX,cellY,'.','color',cColors{ss});
        end
    end
end

nPixels = round(50/micronsPerPixel);

scaleBarS = xrange(2)- 20 - nPixels; scaleBarE = scaleBarS + nPixels;
scaleBarY = yrange(1) + 30;
plot([scaleBarS,scaleBarE],[scaleBarY,scaleBarY],'c','linewidth',3);

% text(scaleBarS,scaleBarY-20,sprintf('%.0f um',(scaleBarE - scaleBarS)*micronsPerPixel),'color','r','FontSize',14,'FontWeight','Bold');
if iscell(selCells)
%     title(tcColors');
else
%     title(tei.recordingFolder)
end
xlim(xrange); ylim(yrange);
disp(tcColors);