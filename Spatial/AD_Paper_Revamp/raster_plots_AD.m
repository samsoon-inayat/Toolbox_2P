an = 4; cn = 1
plotRasters_simplest(Rs_C{an,cn})

%%
magfac = mData.magfac;
ff = makeFigureRowsCols(108,[10 3 5 1.25],'RowsCols',[1 5],'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.07 0.26],'widthHeightAdjustment',[10 -410]);
MY = 10; ysp = 0.5; mY = 0; titletxt = ''; ylabeltxt = {'Trial #'};
stp = 0.3*magfac; widths = ([1.3 1.3 1.3 1.3 1.3 0.5 0.5 0.5]-0.51)*magfac; gap = 0.16*magfac;
adjust_axes(ff,[mY MY],stp,widths,gap,{''});

an = 4; cn = 1;
R = Rs_C{an,cn};
cellN = [322 233 306 313 255];
cbar_p_shift = [-0.011 0.09 -0.03 -0.3];
for cc = 1:5
    c = cellN(cc);
    thisRaster = R.sp_rasters(:,:,c);
    ax = ff.h_axes(1,cc);
    axes(ax); xlabel(ax,'Distance (cm)');
    m = min(thisRaster(:));
    M = max(thisRaster(:));
    xdata = [0 75 150]; ydata = [1 10];
    imagesc(xdata,ydata,thisRaster,[m M]);
    set(gca,'Ydir','normal');
%     mSig = nanmean(thisRaster);
%     m = round(min(mSig),1);
%     M = round(max(mSig),1);
%     hold on;
%     plot(10*mSig/M,'w');
    box off;
    if cc > 1
            set(gca,'ytick',[]);
    else
        ylabel('Trial #');
    end
    
    xlabel('Distance (cm)'); 
    format_axes_b(gca)
    colormap jet;
    textstr = sprintf('%d',c);
    ht = set_axes_top_text_no_line(ff.hf,gca,textstr,[0 -0.07 0 0]); set(ht,'Fontsize',7);
    mM = min(thisRaster(:)); MM = max(thisRaster(:)); [hc,hca] = putColorBar(gca,cbar_p_shift,[mM MM],5,'eastoutside',[0.09 0.11 0.09 0.16]);
end

save_pdf(ff.hf,mData.pdf_folder,'rasters.pdf',600);

