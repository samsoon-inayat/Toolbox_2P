function general_imagesc(fig_props,graphs)

fignumi = fig_props{1};
figpos = fig_props{2};
% 
% 
% ff = makeFigureRowsCols(fignum,[1 1.5 2 2],'RowsCols',[1 numel(popVs)],...
%     'spaceRowsCols',[0.01 0.011],'rightUpShifts',[0.01 0.12],'widthHeightAdjustment',...
%     [-15 -205]);
% gg = 1;
% set(gcf,'color','w');
% set(gcf,'Position',figpos);
figure(fignumi);clf
for rr = 1:size(graphs,1)
    for cc = 1:size(graphs,2)
        subplot(size(graphs,2),size(graphs,1),sub2ind(size(graphs),rr,cc));
        imagesc(graphs{rr,cc});colorbar;
    end
end
set(gcf,'Units','Inches');set(gcf,'Position',figpos)
% for cc = 1:length(popVs)
%     axes(ff.h_axes(1,cc));
%     thisImage = popVs{cc}.popV;
%     mdata = rasters{cc}.mdata;
%     imagesc(thisImage);hold on;
%     if ~isempty(mdata)
%         ys = size(thisImage,1);
%         plot([mdata.cis(1) mdata.cis(1)],[0 ys],'linewidth',1.5,'color','r');
%         plot([mdata.cis(2) mdata.cis(2)],[0 ys],'linewidth',1.5,'color','m');
%     end
%     box off;
%     set(gca,'Ydir','Normal','linewidth',1,'FontSize',6,'FontWeight','Bold');
% end