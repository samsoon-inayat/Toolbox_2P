function raster_plots

%%
%%
an = 4; cn = 3;
% respC = cell_list_op(FD_conj{2},dzMI_FD.resp_complex,'and'); %all_responsive_cells{an,cn}
respC = FD_Dis_comp{2}; respC = FD_conj{2};
figure(1000);clf;subplot 141;imagesc(RsTt{an,cn}.speed); set(gca,'Ydir','normal'); subplot 142;imagesc(RsDt{an,cn}.speed);set(gca,'Ydir','normal'); subplot 143;imagesc(RsTi{an,cn}.speed); set(gca,'Ydir','normal'); subplot 144;imagesc(RsDi{an,cn}.speed);set(gca,'Ydir','normal');
plotRasters_dis_dur({RsTt{an,cn},RsDt{an,cn},RsTi{an,cn},RsDi{an,cn}},find(respC{an,cn}));

%%
an = 4; cn = 3;
Rs = {RsTt{an,cn},RsTi{an,cn}};
all_cellN = [35 32 41 21 86 115 199];
cellN = all_cellN(1);
ff = makeFigureRowsCols(2020,[0.5 0.5 4 1],'RowsCols',[1 2],...
    'spaceRowsCols',[0.15 0.08],'rightUpShifts',[0.09 0.25],'widthHeightAdjustment',...
    [-100 -475]);
set(gcf,'color','w'); set(gcf,'Position',[10 4 3 1]);
ff = sample_rasters(Rs,cellN,ff);
axes(ff.h_axes(1,1));ylabel('Trial #');
save_pdf(ff.hf,mData.pdf_folder,sprintf('air_rastersT'),600);

Rs = {RsDt{an,cn},RsDi{an,cn}};
ff = makeFigureRowsCols(2021,[0.5 0.5 4 1],'RowsCols',[1 2],...
    'spaceRowsCols',[0.15 0.08],'rightUpShifts',[0.09 0.25],'widthHeightAdjustment',...
    [-100 -475]);
set(gcf,'color','w'); set(gcf,'Position',[10 7 3 1]);
ff = sample_rasters(Rs,cellN,ff);
axes(ff.h_axes(1,1));ylabel('Trial #');
save_pdf(ff.hf,mData.pdf_folder,sprintf('air_rastersD'),600);


