function raster_plots

%%
%%
an = 4; cn = 2;
% respC = cell_list_op(FD_conj{2},dzMI_FD.resp_complex,'and'); %all_responsive_cells{an,cn}
respC = FD_Dis_comp{3};  respC = FD_conj{3};
respC = repmat(cell_list_op(dis_cells_TC,dur_cells_IC,'and'),1,3);
respC = repmat(dis_cells_I(:,cn),1,3);
p_respC = find_percent(respC);

figure(1000);clf;subplot 141;imagesc(RsTt{an,cn}.speed); set(gca,'Ydir','normal'); subplot 142;imagesc(RsDt{an,cn}.speed);set(gca,'Ydir','normal'); subplot 143;imagesc(RsTi{an,cn}.speed); set(gca,'Ydir','normal'); subplot 144;imagesc(RsDi{an,cn}.speed);set(gca,'Ydir','normal');
cNs = find(respC{an,cn});
% cNs = find(speedResp{an,cn});
% cNs = [278 118 188 97 35 329 115 85 21 209 132 238 251 149];
% cNs = [278 97 329 209 132 251 149];
tprops = get_props_Rs(RsTi(an,cn),50);
tprops = get_props_Rs(RsDt(an,cn),50);
centers = tprops.centers{1}(cNs,1); [ic,vc] = sort(centers);
cNs = cNs(vc);
cNs = find(respC{an,cn});
% cNs = cNs(find(mean_dzMI{1,10}'==0))'
plotRasters_dis_dur({RsTt{an,cn},RsDt{an,cn},RsTi{an,cn},RsDi{an,cn}},cNs);

%%
an = 4; cn = 2;
tRsTt = RsTtMC{an,cn};tRsTi = RsTiMC{an,cn};
tRsDt = RsDtMC{an,cn};tRsDi = RsDiMC{an,cn};
all_cellN = [207 123 148 329 335 108 264];
all_cellN = [278 97 329 209 132 251 149];
all_cellN = [21 32 145];
all_cellN = [6 36 303 140];
all_cellN = 1;
for ii = 1:length(all_cellN)
    cellN = all_cellN(ii);
    ff = makeFigureRowsCols(2020,[0.5 0.5 4 1],'RowsCols',[1 4],...
        'spaceRowsCols',[0.15 0.05],'rightUpShifts',[0.05 0.25],'widthHeightAdjustment',...
        [-70 -800]);
    set(gcf,'color','w'); set(gcf,'Position',[10 7 6.9 1.15]); sh34 = 0.06;
    ax = ff.h_axes(1,3); changePosition(ax,[sh34 0 0 0]); ax = ff.h_axes(1,4); changePosition(ax,[sh34 0 0 0]);
    
    ax = ff.h_axes(1,1); [hra,hca] = plot_raster(tRsTt,cellN,ax); xlabel(ax,'Time (s)'); 
    th = ylabel(ax,'FR (A.U.)'); changePosition(th,[-0.5,0,0]);
    th = ylabel(hra,'Trial #'); changePosition(th,[-2.25,0,0]);
    ax = ff.h_axes(1,2); [hra,hca] = plot_raster(tRsDt,cellN,ax); xlabel(ax,'Distance (cm)');
    ax = ff.h_axes(1,3); [hra,hca] = plot_raster(tRsTi,cellN,ax); xlabel(ax,'Time (s)');
    ax = ff.h_axes(1,4); [hra,hca] = plot_raster(tRsDi,cellN,ax); xlabel(ax,'Distance (cm)');
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


