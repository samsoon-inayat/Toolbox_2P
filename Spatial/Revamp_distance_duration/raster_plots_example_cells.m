an = 2; cn = 15;
plotRasters_simplest(Rs_C{an,cn},find(sel_pop_C{an,cn}))
%% rasters
magfac = mData.magfac;
ff = makeFigureRowsCols(108,[10 3 6.9 4.5],'RowsCols',[3 6],'spaceRowsCols',[0.19 -0.02],'rightUpShifts',[0.02 0.1],'widthHeightAdjustment',[10 -200]);
MY = 10; ysp = 0.5; mY = 0; titletxt = ''; ylabeltxt = {'Trial #'};
stp = 0.6*magfac; widths = (ones(1,18)*0.85)*magfac; gap = 0.2*magfac;
adjust_axes(ff,[mY MY],stp,widths,gap,{''});

zvalth = [1.65 1.65 1.65];
si_cn_ap_H = si_cn_ap(:,[1 3 5 2 4 6]);


%  time
an = 4; 
condNum = 6;
mod_cells = propsTDM.newMI.cells_time{4,condNum};
zvals = propsTDM.newMI.cells_time_zvals{4,condNum};
zvalsPC = propsTDM.newPC.cells_time_zvals{4,condNum};
cellN = find(mod_cells & zvals(:,1) > zvalth(1) & abs(zvalsPC(:,1)) > zvalth(1)); cellN = setdiff(cellN,[58,57,192,7,22,77]);
zvalsV = zvals(cellN,1); zvalsVPC = zvalsPC(cellN,1);
si = [Ars_i_T Ars_i_T Ars_i_T Ars_i_T Ars_i_T Ars_i_T];
Rs_C = o.Rs(:,si);mRs_C = o.mR(:,si);
Rs = Rs_C(an,:);
cbar_p_shift = [0.01 0.09 -0.05 -0.3];
R = Rs{condNum};
MVT = met_valsT{1}{an,si_cn_ap_H(1,condNum),si_cn_ap_H(2,condNum)};
MVD = met_valsD{1}{an,si_cn_ap_H(1,condNum),si_cn_ap_H(2,condNum)};
for cc = 1:6
    c = cellN(cc);
    thisRaster = R.sp_rasters(:,:,c);
    zval_here = [MVT.PC(c,2) MVT.MI(c,2)];
    zval_here = [zvalsVPC(cc) zvalsV(cc)];
    ax = ff.h_axes(1,cc);
    axes(ax); xlabel(ax,'Time (s)');
    m = min(thisRaster(:));
    M = max(thisRaster(:));
    sz = size(thisRaster,2);
    all_sz(cc) = sz;
    bw = R.bin_width;
    xdata = [0 round(sz*bw/2) sz*bw]; ydata = [1 10];
    try
    imagesc(xdata,ydata,thisRaster,[m M]);
    catch
      continue;
    end
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
    try
    [hc,hca] = putColorBar(gca,cbar_p_shift,[mM MM],5,'eastoutside',[0.1 0.11 0.09 0.16]);
    catch
      continue;
    end

    %******* make new axes and plot mean and gaussian fitting
    pos = get(ax,'Position'); ha = axes; set(ha,'units','inches');set(ha,'Position',pos);
    changePosition(ha,[0 pos(4)+0.05 0 -0.35]);

    xs = linspace(0,sz*bw,sz);
    mSig = nanmean(thisRaster);
    plot(xs,mSig,'b');hold on;

    fitplot = gauss_fit(1:length(xs),R.gauss_fit_on_mean.coefficients_Rs_mean(c,1:3),R.gauss_fit_on_mean.gauss1Formula);
    % plot(xs,fitplot,'linewidth',0.5,'color','m');
    box off;
    try
    % ylim([0 round(max(mSig)+(max(mSig)/10),1)])
    catch
      continue;
    end
    set(ha,'xtick',[],'ytick',round(max(mSig),1));
    format_axes_b(gca);
    textstr = sprintf('Cell %d (%.1f, %.1f)',c,zvalsV(cc),R.gauss_fit_on_mean.coefficients_Rs_mean(c,4));
    % textstr = sprintf('%.1f',zvalsV(cc));
    textstr = sprintf('Cell-%d, (%.1f, %.1f)',c,zval_here);
    ht = set_axes_top_text_no_line(ff.hf,gca,textstr,[0.01 0.02 0 0]); set(ht,'Fontsize',6,'FontWeight','normal');
    if cc == 1
        ylabel('FR (AU)');
    end
end

% dist
condNum = 5;
mod_cells = propsTDM.newMI.cells_dist{4,condNum};
zvals = propsTDM.newMI.cells_dist_zvals{4,condNum};
zvalsPC = propsTDM.newPC.cells_dist_zvals{4,condNum};
cellN = find(mod_cells & zvals(:,2) > zvalth(2) & abs(zvalsPC(:,2)) > zvalth(2)); cellN = setdiff(cellN,[312,201,122]);
zvalsV = zvals(cellN,2); zvalsVPC = zvalsPC(cellN,2);
si = [Ar_t_D Ar_i_D ArL_t_D ArL_i_D Ars_t_D Ars_i_D];
Rs_C = o.Rs(:,si);mRs_C = o.mR(:,si);
cbar_p_shift = [0.01 0.09 -0.05 -0.3];
Rs = Rs_C(an,:);
cbar_p_shift = [0.01 0.09 -0.05 -0.3];
R = Rs{condNum};
for cc = 1:6
    c = cellN(cc);
    thisRaster = R.sp_rasters(:,:,c);
    zval_here = [MVD.PC(c,2) MVD.MI(c,2)];
    zval_here = [zvalsVPC(cc) zvalsV(cc)];
    ax = ff.h_axes(2,cc);
    axes(ax); xlabel(ax,'Distance (cm)');
    m = min(thisRaster(:));
    M = max(thisRaster(:));
    sz = size(thisRaster,2);
    all_sz(cc) = sz;
    bw = R.bin_width;
    xdata = [0 round(sz*bw/2) sz*bw]; ydata = [1 10];
    try
      imagesc(xdata,ydata,thisRaster,[m M]);
    catch
      continue;
    end
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
    try
    [hc,hca] = putColorBar(gca,cbar_p_shift,[mM MM],5,'eastoutside',[0.1 0.11 0.09 0.16]);
    catch
      continue;
    end

    %******* make new axes and plot mean and gaussian fitting
    pos = get(ax,'Position'); ha = axes; set(ha,'units','inches');set(ha,'Position',pos);
    changePosition(ha,[0 pos(4)+0.05 0 -0.35]);

    xs = linspace(0,sz*bw,sz);
    mSig = nanmean(thisRaster);
    plot(xs,mSig,'b');hold on;
    fitplot = gauss_fit(1:length(xs),R.gauss_fit_on_mean.coefficients_Rs_mean(c,1:3),R.gauss_fit_on_mean.gauss1Formula);
    % plot(xs,fitplot,'linewidth',0.5,'color','m');
    box off;
    try
    % ylim([0 round(max(mSig)+(max(mSig)/10),1)])
    catch
      continue;
    end
    set(ha,'xtick',[],'ytick',round(max(mSig),1));
    format_axes_b(gca);
%     textstr = sprintf('Cell %d (%.1f, %.1f)',c,R.info_metrics.ShannonMI_Zsh(c),R.gauss_fit_on_mean.coefficients_Rs_mean(c,4));
    % textstr = sprintf('%.1f',zvalsV(cc));
    textstr = sprintf('Cell-%d, (%.1f, %.1f)',c,zval_here);
    ht = set_axes_top_text_no_line(ff.hf,gca,textstr,[0.01 0.02 0 0]); set(ht,'Fontsize',6,'FontWeight','normal');
    if cc == 1
        ylabel('FR (AU)');
    end
end

%  time
an = 4; 
condNum = 6;
mod_cells = propsTDM.newPC.cells_speed{4,condNum};
zvals = propsTDM.newMI.cells_speed_zvals{4,condNum};
zvalsPC = propsTDM.newPC.cells_speed_zvals{4,condNum};
cellN = find(mod_cells & zvals(:,3) > zvalth(3) & abs(zvalsPC(:,3)) > zvalth(3)); cellN = setdiff(cellN,[3,35,22,33,62,68,65]);
zvalsV = zvals(cellN,3); zvalsVPC = zvalsPC(cellN,3);
si = [Ars_i_T Ars_i_T Ars_i_T Ars_i_T Ars_i_T Ars_i_T];
Rs_C = o.Rs(:,si);mRs_C = o.mR(:,si);
Rs = Rs_C(an,:);
cbar_p_shift = [0.01 0.09 -0.05 -0.3];
R = Rs{condNum};

d.bcs = speedRs{an}.bin_centers;
d.FR = speedRs{an}.FR_vs_speed;

for cc = 1:6
    c = cellN(cc);
    thisRaster = R.sp_rasters(:,:,c);
    zval_here = [MVT.PC(c,2) MVT.MI(c,2)];
    zval_here = [zvalsVPC(cc) zvalsV(cc)];
    ax = ff.h_axes(3,cc);
    axes(ax); 
    % plot(d.bcs,d.FR(c,:));
    % continue;
    xlabel(ax,'Time (s)');
    m = min(thisRaster(:));
    M = max(thisRaster(:));
    sz = size(thisRaster,2);
    all_sz(cc) = sz;
    bw = R.bin_width;
    xdata = [0 round(sz*bw/2) sz*bw]; ydata = [1 10];
    try
    imagesc(xdata,ydata,thisRaster,[m M]);
    catch
      continue;
    end
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
    try
    [hc,hca] = putColorBar(gca,cbar_p_shift,[mM MM],5,'eastoutside',[0.1 0.11 0.09 0.16]);
    catch
      continue;
    end

    %******* make new axes and plot mean and gaussian fitting
    pos = get(ax,'Position'); ha = axes; set(ha,'units','inches');set(ha,'Position',pos);
    changePosition(ha,[0 pos(4)+0.05 0 -0.35]);

    xs = linspace(0,sz*bw,sz);
    mSig = nanmean(thisRaster);
    mSpd = nanmean(R.speed);
    axyy = plotyy(xs,mSig,xs,mSpd);hold on;
    set(axyy(2),'ylim',[0 40])
    if cc < 6
        set(axyy(2),'ytick',[]);
    else
        set(axyy(2),'xtick',[],'ytick',40);
    end

    fitplot = gauss_fit(1:length(xs),R.gauss_fit_on_mean.coefficients_Rs_mean(c,1:3),R.gauss_fit_on_mean.gauss1Formula);
    % plot(xs,fitplot,'linewidth',0.5,'color','m');
    box off;
    try
    % ylim([0 round(max(mSig)+(max(mSig)/10),1)])
    catch
      continue;
    end
    set(axyy(1),'xtick',[],'ytick',round(max(mSig),1));
    
    format_axes_b(axyy(2))
    format_axes_b(axyy(1));
    set(axyy(2),'ycolor','r')
    textstr = sprintf('Cell-%d, (%.1f, %.1f)',c,zval_here);
    ht = set_axes_top_text_no_line(ff.hf,gca,textstr,[0.01 0.02 0 0]); set(ht,'Fontsize',6,'FontWeight','normal');
    if cc == 1
        ylabel('FR (AU)');
    end
end



save_pdf(ff.hf,mData.pdf_folder,'rasters.pdf',600);


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

