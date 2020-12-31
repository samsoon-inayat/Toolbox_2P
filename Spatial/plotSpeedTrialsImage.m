function plotTrialsImage(b,markers1,markers2,fn)
%%
n = 0;
%%
ei = evalin('base','ei10_A');
% ei = ei([1:4 9]);
mData = evalin('base','mData');
ii = 4; cc = 4; selAnimal = ii;
b = ei{ii}.b;
markers1i = ei{ii}.plane{1}.contexts(cc).markers.air_onsets;
markers2i = ei{ii}.plane{1}.contexts(cc).markers.air_offsets;
fn = 101;

speed = b.fSpeed;
% figure(1000);clf;plot(b.speed);
ts = b.ts;

timeBefore = 0;
timeAfter = 15;

markers1 = markers1i - round(1e6 * timeBefore/b.si);
markers2 = markers2i + round(1e6 * timeAfter/b.si);

for ii = 1:length(markers1)
    st = markers1(ii);
    se = markers2(ii);
    sp{ii} = speed(st:se);
    lsp(ii) = length(sp{ii});
    t{ii} = ts(st:se)-ts(st);
    ind(ii) = find((st:se)-markers2i(ii)>0,1,'first');
end

mlsp = min(lsp);
for ii = 1:length(markers1)
    isp(ii,:) = sp{ii}(1:mlsp);
end
it = t{1}(1:mlsp);

% ff = makeFigureWindow__one_axes_only(5,[6 4 5 5],[0.1 0.012 0.85 0.99]);
figure(101);clf;
set(gcf,'units','inches');
set(gcf,'Position',[15 4 1 1]);
set(gcf,'color','w');
% for ii = 1:5
%     aisp(:,:,ii) = ei{ii}.plane{1}.contexts(1).rasters.airD.speed;
% end
% % isp = ei{8}.plane{1}.contexts(1).rasters.airT.speed;
% isp = mean(aisp,3);
% ispsem = std(aisp,[],3)/sqrt(4);
imagesc(isp);

% hc = colorbar;changePosition(hc,[0.19 -0.1 -0.02 0.02]);
% set(hc,'linewidth',0.25,'TickDirection','out');
changePosition(gca,[-0.02 -0.2 -0.25 0.07]);
box off
xtt = 0:5:20;
for ii = 1:length(xtt)
    xt(ii) = find(it>xtt(ii),1,'first');
end
txt = cellstr(num2str(ceil(xtt)'));

% axis equal
set(gca,'Ydir','normal','XTick',xt,'XTickLabel',txt,'FontSize',6,'FontWeight','Bold','TickDir','out','linewidth',1);
h = ylabel('Trials'); changePosition(h,[6000 0 0]);
h = xlabel('Time (sec)');changePosition(h,[0 0.3 0]);
text(length(it)+50000,1,'Speed (cm/sec)','rotation',90,'FontSize',5);
% n = 0;
hold on;
for ii = 1:length(markers1)
    plot([ind(ii) ind(ii)],[ii-0.5 ii+0.5],'r','linewidth',0.75)
end
changePosition(gca,[0.07 0.08 0 0]);
% [x1 y1] = ds2nfu(-8000,13.5);
% [x2 y2] = ds2nfu(0,11.75);
% annotation('textarrow',[x1 x2],[y1 y2],'String','Air onset','HeadLength',3,'HeadWidth',3,'FontSize',5)
% [x1 y1] = ds2nfu(ind(ii)+6000,13);
% [x2 y2] = ds2nfu(ind(ii),11.6);
% annotation('textarrow',[x1 x2],[y1 y2],'String','Air offset','HeadLength',3,'HeadWidth',3,'FontSize',5,'color','r')
hc = putColorBar(gca,[-0.01 0.09 0 -0.2],[min(isp(:)) max(isp(:))],5,'text_margins',[0.1 0.15]);
colormap parula;

save_pdf(gcf,mData.pdf_folder,sprintf('speedTrialsImage_%d.pdf',cc),600);
return;

