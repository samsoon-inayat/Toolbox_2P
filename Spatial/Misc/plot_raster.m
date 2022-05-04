function [ha,hc] = plot_raster(R,cn,ax,sp)
if ~exist('sp','var')
    sp = 0;
end

A = R;
axes(ax);
if sp
    thisRaster = A.speed;
else
    thisRaster = A.sp_rasters(:,:,cn);
end
    
mSig = nanmean(thisRaster);
m = round(min(mSig),1);
M = round(max(mSig),1);
xs = A.xs;
xs(length(xs)+1) = (xs(2)-xs(1))+xs(end);
fitplot = gauss_fit(1:length(xs),A.gauss_fit_on_mean.coefficients_Rs_mean(cn,1:3),A.gauss_fit_on_mean.gauss1Formula);
box off;
plot(xs,mSig,'b','linewidth',1);hold on;
if ~sp
    plot(xs,fitplot,'linewidth',0.5,'color','m');
end
box off;
cols = size(thisRaster,2);
colsHalf = ceil(cols/2);
xlim([xs(1) xs(end)]);
format_axes(ax);
%% raster
pos = get(ax,'Position');
ha = axes('Position',pos,'Visible','on');
changePosition(ha,[0 0.3 0 0.15]);
m = min(thisRaster(:));
M = max(thisRaster(:));
imagesc(thisRaster,[m M]);
set(gca,'Ydir','normal');
% axis off;
box off;
set(gca,'xtick',[],'ytick',[1 size(thisRaster,1)],'TickDir','out');
format_axes(gca)
ha.XAxis.Visible = 'Off';
if ~sp
text(1,size(thisRaster,1)+2,sprintf('zMI = %.2f ',R.info_metrics.ShannonMI_Zsh(cn)),'FontSize',5,'color','k');
end
hold on;
if ~isempty(strfind(R.context_info,'airT')) || ~isempty(strfind(R.context_info,'airID'))
    for bb = 1:length(R.lastBin)
        xvalbin = R.lastBin(bb);
        plot([xvalbin xvalbin]+0.15,[bb-0.75 bb+0.75],'r','linewidth',0.25);
    end
end
colormap parula;
%% raster colorbar
hca = ha;
MS = sprintf('%.1f',M);
mS = sprintf('%.1f',abs(m));
hc = putColorBar(hca,[0.1 0.112 -0.125 -0.2],{mS,MS},6,'northoutside',[0.4 1 0.05 1]);
format_axes(ax);




