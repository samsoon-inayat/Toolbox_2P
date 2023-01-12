function raster_plots


%%
an = 1; cn = 1;

asi = [Ab_On Abs_Off Ar_On Ar_Off];
    si = asi(cn);
%     si = [Lb];
    Rs = o.RsMC(:,asi);
    
tRsTt = Rs{an,1};tRsTi = Rs{an,3};
tRsDt = Rs{an,2};tRsDi = Rs{an,4};
all_cellN = 1;
for ii = 1:length(all_cellN)
    cellN = all_cellN(ii);
    ff = makeFigureRowsCols(2020,[0.5 0.5 4 1],'RowsCols',[1 4],...
        'spaceRowsCols',[0.15 0.05],'rightUpShifts',[0.05 0.25],'widthHeightAdjustment',...
        [-70 -800]);
    set(gcf,'color','w'); set(gcf,'Position',[10 7 6.9 1.15]); sh34 = 0.06;
    ax = ff.h_axes(1,3); changePosition(ax,[sh34 0 0 0]); ax = ff.h_axes(1,4); changePosition(ax,[sh34 0 0 0]);
    
    ax = ff.h_axes(1,1); [hra,hca] = plot_raster(tRsTt,cellN,ax); xlabel(ax,'Time (s)'); 
    th = ylabel(ax,'MC (pix)'); changePosition(th,[-0.5,0,0]);
    th = ylabel(hra,'Trial #'); changePosition(th,[-2.25,0,0]);
    ax = ff.h_axes(1,2); [hra,hca] = plot_raster(tRsDt,cellN,ax); xlabel(ax,'Time (s)');
    ax = ff.h_axes(1,3); [hra,hca] = plot_raster(tRsTi,cellN,ax); xlabel(ax,'Time (s)');
    ax = ff.h_axes(1,4); [hra,hca] = plot_raster(tRsDi,cellN,ax); xlabel(ax,'Time (s)');
    save_pdf(ff.hf,mData.pdf_folder,sprintf('air_rasters_%d',cellN),600);
end

%%
for ii = 1%:length(all_cellN)
    cellN = 1;%all_cellN(ii);
    ff = makeFigureRowsCols(2020,[0.5 0.5 4 1],'RowsCols',[1 4],...
        'spaceRowsCols',[0.15 0.05],'rightUpShifts',[0.05 0.25],'widthHeightAdjustment',...
        [-70 -800]);
    set(gcf,'color','w'); set(gcf,'Position',[10 7 6.9 1.15]); sh34 = 0.06;
    ax = ff.h_axes(1,3); changePosition(ax,[sh34 0 0 0]); ax = ff.h_axes(1,4); changePosition(ax,[sh34 0 0 0]);
    
    ax = ff.h_axes(1,1); [hra,hca] = plot_raster(tRsTt,cellN,ax,1); xlabel(ax,'Time (s)'); 
    th = ylabel(ax,'(cm/sec)'); changePosition(th,[-0.5,0,0]);
    th = ylabel(hra,'Trial #'); changePosition(th,[-2.25,0,0]);
    ax = ff.h_axes(1,2); [hra,hca] = plot_raster(tRsDt,cellN,ax,1); xlabel(ax,'Distance (cm)');
    ax = ff.h_axes(1,3); [hra,hca] = plot_raster(tRsTi,cellN,ax,1); xlabel(ax,'Time (s)');
    ax = ff.h_axes(1,4); [hra,hca] = plot_raster(tRsDi,cellN,ax,1); xlabel(ax,'Distance (cm)');
    save_pdf(ff.hf,mData.pdf_folder,sprintf('speed_%d',1),600);
end
%%
axes(ff.h_axes(1,1));ylabel('Trial #');

Rs = {RsDt{an,cn},RsDi{an,cn}};
ff = makeFigureRowsCols(2021,[0.5 0.5 4 1],'RowsCols',[1 2],...
    'spaceRowsCols',[0.15 0.06],'rightUpShifts',[0.09 0.25],'widthHeightAdjustment',...
    [-100 -475]);
set(gcf,'color','w'); set(gcf,'Position',[10 7 3 1]);
ff = sample_rasters(Rs,cellN,ff);
axes(ff.h_axes(1,1));ylabel('Trial #');
save_pdf(ff.hf,mData.pdf_folder,sprintf('air_rastersD'),600);


