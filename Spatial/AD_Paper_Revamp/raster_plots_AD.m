an = 4; cn = 3;
plotRasters_simplest(Rs_C{an,cn},find(sel_pop_C{an,cn}))
%%
an = 2; cn = 3;
plotRasters_simplest(Rs_A{an,cn},find(sel_pop_A{an,cn}))

%%
while 1
    magfac = mData.magfac;
    ff = makeFigureRowsCols(108,[10 3 3.5 1.25],'RowsCols',[1 4],'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.02 0.20],'widthHeightAdjustment',[10 -550]);
    MY = 10; ysp = 0.5; mY = 0; titletxt = ''; ylabeltxt = {'Trial #'};
    stp = 0.27*magfac; widths = ([1.3 1.3 1.3 1.3 1.3 0.5 0.5 0.5]-0.71)*magfac; gap = 0.24*magfac;
    adjust_axes(ff,[mY MY],stp,widths,gap,{''});

    an = 4; cn = 3;R = Rs_C{an,cn};cellN = [68 33 255 88 181];cellN = [212 288 275 299];
    an = 2; cn = 3;R = Rs_A{an,cn};cellN = [222 199 5 73];
    cbar_p_shift = [0.01 0.09 -0.05 -0.3];
    for cc = 1:4
        c = cellN(cc);
        thisRaster = R.sp_rasters(:,:,c);
        ax = ff.h_axes(1,cc);
        axes(ax); xlabel(ax,'Distance (cm)');
        m = min(thisRaster(:));
        M = max(thisRaster(:));
        xdata = [0 75 150]; ydata = [1 10];
        imagesc(xdata,ydata,thisRaster,[m M]);
        hold on;
        set(gca,'Ydir','normal');
        if cc > 1
                set(gca,'ytick',[]);
        else
            ylabel('Trial #');
        end
        xlabel('Distance (cm)'); 
        format_axes_b(gca)
        colormap parula;
        box off;

        %****** color bar
        mM = round(min(thisRaster(:)),1); MM = round(max(thisRaster(:)),1);
        [hc,hca] = putColorBar(gca,cbar_p_shift,[mM MM],5,'eastoutside',[0.1 0.11 0.09 0.16]);

        %******* make new axes and plot mean and gaussian fitting
        pos = get(ax,'Position'); ha = axes; set(ha,'units','inches');set(ha,'Position',pos);
        changePosition(ha,[0 pos(4)+0.05 0 -0.35]);

        xs = linspace(0,150,50);
        mSig = nanmean(thisRaster);
        plot(xs,mSig,'b');hold on;
        fitplot = gauss_fit(1:length(xs),R.gauss_fit_on_mean.coefficients_Rs_mean(c,1:3),R.gauss_fit_on_mean.gauss1Formula);
        plot(xs,fitplot,'linewidth',0.5,'color','m');
        box off;
        ylim([0 round(max(mSig)+(max(mSig)/10),1)])
        set(ha,'xtick',[],'ytick',round(max(mSig),1));
        format_axes_b(gca);
        textstr = sprintf('Cell %d (%.1f, %.1f)',c,R.info_metrics.ShannonMI_Zsh(c),R.gauss_fit_on_mean.coefficients_Rs_mean(c,4));
        textstr = sprintf('Cell %d (%.1f)',c,R.info_metrics.ShannonMI_Zsh(c));
        ht = set_axes_top_text_no_line(ff.hf,gca,textstr,[-0.01 -0.05 0 0]); set(ht,'Fontsize',6,'FontWeight','Bold');
        if cc == 1
            ylabel('FR (AU)');
        end
    end

    save_pdf(ff.hf,mData.pdf_folder,'rasters.pdf',600);
    %%
break;
end

%% Time rasters
while 1
    magfac = mData.magfac;
    ff = makeFigureRowsCols(108,[10 3 3.5 1.25],'RowsCols',[1 4],'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.02 0.20],'widthHeightAdjustment',[10 -550]);
    MY = 10; ysp = 0.5; mY = 0; titletxt = ''; ylabeltxt = {'Trial #'};
    stp = 0.27*magfac; widths = ([1.3 1.3 1.3 1.3 1.3 0.5 0.5 0.5]-0.71)*magfac; gap = 0.24*magfac;
    adjust_axes(ff,[mY MY],stp,widths,gap,{''});

    an = 4; cn = 3;R = Rs_C{an,cn};cellN = [155 37 246 92]; %21 51
%     an = 2; cn = 3;R = Rs_A{an,cn};cellN = [15 191 72 155 95 64 ];
    cbar_p_shift = [0.01 0.09 -0.05 -0.3];
    for cc = 1:4
        c = cellN(cc);
        thisRaster = R.sp_rasters(:,1:50,c);

        ax = ff.h_axes(1,cc);
        axes(ax); xlabel(ax,'Time (s)');
        m = min(thisRaster(:));
        M = max(thisRaster(:));
        xdata = [0 7.5 15]; ydata = [1 10];
        imagesc(xdata,ydata,thisRaster,[m M]);
        hold on;
        set(gca,'Ydir','normal');
        if cc > 1
                set(gca,'ytick',[]);
        else
            ylabel('Trial #');
        end
        xlabel('Time (s)'); 
        format_axes_b(gca)
        colormap parula;
        box off;

        %****** color bar
        mM = round(min(thisRaster(:)),1); MM = round(max(thisRaster(:)),1);
        [hc,hca] = putColorBar(gca,cbar_p_shift,[mM MM],5,'eastoutside',[0.1 0.11 0.09 0.16]);

        %******* make new axes and plot mean and gaussian fitting
        pos = get(ax,'Position'); ha = axes; set(ha,'units','inches');set(ha,'Position',pos);
        changePosition(ha,[0 pos(4)+0.05 0 -0.35]);

        xs = linspace(0,150,50);
        mSig = nanmean(thisRaster);
        plot(xs,mSig,'b');hold on;
        fitplot = gauss_fit(1:length(xs),R.gauss_fit_on_mean.coefficients_Rs_mean(c,1:3),R.gauss_fit_on_mean.gauss1Formula);
        plot(xs,fitplot,'linewidth',0.5,'color','m');
        box off;
        ylim([0 round(max(mSig)+(max(mSig)/10),1)])
        set(ha,'xtick',[],'ytick',round(max(mSig),1));
        format_axes_b(gca);
        textstr = sprintf('Cell %d (%.1f, %.1f)',c,R.info_metrics.ShannonMI_Zsh(c),R.gauss_fit_on_mean.coefficients_Rs_mean(c,4));
        textstr = sprintf('Cell %d (%.1f)',c,R.info_metrics.ShannonMI_Zsh(c));
        ht = set_axes_top_text_no_line(ff.hf,gca,textstr,[-0.01 -0.05 0 0]); set(ht,'Fontsize',6,'FontWeight','Bold');
        if cc == 1
            ylabel('FR (AU)');
        end
    end

    save_pdf(ff.hf,mData.pdf_folder,'rasters.pdf',600);
    %%
break;
end

%%
magfac = mData.magfac;
ff = makeFigureRowsCols(108,[10 3 5 1.5],'RowsCols',[1 5],'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.02 0.16],'widthHeightAdjustment',[10 -500]);
MY = 10; ysp = 0.5; mY = 0; titletxt = ''; ylabeltxt = {'Trial #'};
stp = 0.3*magfac; widths = ([1.3 1.3 1.3 1.3 1.3 0.5 0.5 0.5]-0.61)*magfac; gap = 0.25*magfac;
adjust_axes(ff,[mY MY],stp,widths,gap,{''});

an = 4; cn = 1;
R = Rs_C{an,cn};
cellN = [68 33 255 88 181];
cbar_p_shift = [0.01 0.09 -0.05 -0.3];
for cc = 1:5
    c = cellN(cc);
    thisRaster = R.sp_rasters(:,:,c);
    ax = ff.h_axes(1,cc);
    axes(ax); xlabel(ax,'Distance (cm)');
    m = min(thisRaster(:));
    M = max(thisRaster(:));
    xdata = [0 75 150]; ydata = [1 10];
    imagesc(xdata,ydata,thisRaster,[m M]);
    hold on;
    set(gca,'Ydir','normal');
    if cc > 1
            set(gca,'ytick',[]);
    else
        ylabel('Trial #');
    end
    xlabel('Distance (cm)'); 
    format_axes_b(gca)
    colormap parula;
    box off;
    
    %****** color bar
    mM = round(min(thisRaster(:)),1); MM = round(max(thisRaster(:)),1);
    [hc,hca] = putColorBar(gca,cbar_p_shift,[mM MM],5,'eastoutside',[0.1 0.11 0.09 0.16]);
    
    %******* make new axes and plot mean and gaussian fitting
    pos = get(ax,'Position'); ha = axes; set(ha,'units','inches');set(ha,'Position',pos);
    changePosition(ha,[0 pos(4)+0.05 0 -0.5]);
    
    xs = linspace(0,150,50);
    mSig = nanmean(thisRaster);
    plot(xs,mSig,'b');hold on;
    fitplot = gauss_fit(1:length(xs),R.gauss_fit_on_mean.coefficients_Rs_mean(c,1:3),R.gauss_fit_on_mean.gauss1Formula);
    plot(xs,fitplot,'linewidth',0.5,'color','m');
    box off;
    ylim([0 round(max(mSig),1)])
    set(ha,'xtick',[],'ytick',round(max(mSig),1));
    format_axes_b(gca);
    textstr = sprintf('Cell %d (%.1f, %.1f)',c,R.info_metrics.ShannonMI_Zsh(c),R.gauss_fit_on_mean.coefficients_Rs_mean(c,4));
    ht = set_axes_top_text_no_line(ff.hf,gca,textstr,[-0.01 -0.05 0 0]); set(ht,'Fontsize',6,'FontWeight','Bold');
    if cc == 1
        ylabel('FR (AU)');
    end
end

save_pdf(ff.hf,mData.pdf_folder,'rasters.pdf',600);

