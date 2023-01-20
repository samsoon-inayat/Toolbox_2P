an = 4;
si = [Lb Ab_On Ab_Off Ar_On Ar_Off Lbs Abs_On Abs_Off];
Rs = o.RsMC(:,si);
% props1 = get_props_Rs(Rs,ntrials);
allR = Rs(an,:);
event_types = rasterNamesTxt(si)

%% Time rasters

    magfac = mData.magfac;
    ff = makeFigureRowsCols(108,[10 3 6.9 1.3],'RowsCols',[1 8],'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.02 0.20],'widthHeightAdjustment',[10 -580]);
    MY = 10; ysp = 0.5; mY = 0; titletxt = ''; ylabeltxt = {'Trial #'};
    stp = 0.27*magfac; widths = ([ones(1,12)*1]-0.41)*magfac; gap = 0.24*magfac;
    adjust_axes(ff,[mY MY],stp,widths,gap,{''});

    cellN = [1]; %21 51
%     an = 2; cn = 3;R = Rs_A{an,cn};cellN = [15 191 72 155 95 64 ];
    cbar_p_shift = [0.01 0.09 -0.05 -0.3];
    for cc = 1:length(si)
        c = 1;
        R = allR{cc};
        thisRaster = R.sp_rasters(:,:,c)*pixelSize(an);

        ax = ff.h_axes(1,cc);
        axes(ax); xlabel(ax,'Time (s)');
        m = min(thisRaster(:));
        M = max(thisRaster(:));
        xdata = [-2 0 2]; ydata = [1 10];
        imagesc(xdata,ydata,thisRaster,[m M]);
        hold on;
        set(gca,'Ydir','normal');
        if cc > 1
                set(gca,'ytick',[]);
        else
            ylabel('Trial #');
        end
        xlabel('Time (s)'); 
        format_axes(gca)
        colormap parula;
        box off;
        ylims = ylim;
        plot([0 0],ylims,'m');
        %****** color bar
        mM = round(min(thisRaster(:)),1); MM = round(max(thisRaster(:)),1);
        [hc,hca] = putColorBar(gca,cbar_p_shift,[mM MM],5,'eastoutside',[0.1 0.11 0.09 0.16]);

        %******* make new axes and plot mean and gaussian fitting
        pos = get(ax,'Position'); ha = axes; set(ha,'units','inches');set(ha,'Position',pos);
        changePosition(ha,[0 pos(4)+0.05 0 -0.35]);

        xs = linspace(-2,2,size(thisRaster,2));
        mSig = nanmean(thisRaster);
        plot(xs,mSig,'b');hold on;
        fitplot = gauss_fit(1:length(xs),R.gauss_fit_on_mean.coefficients_Rs_mean(c,1:3),R.gauss_fit_on_mean.gauss1Formula);
%         plot(xs,fitplot,'linewidth',0.5,'color','m');
        box off;
        ylim([0 round(max(mSig)+(max(mSig)/10),1)])
        ylim([0 8]);
        set(ha,'xtick',[],'ytick',[0 8]);
        format_axes(gca);
        ylims = ylim;
        plot([0 0],ylims,'m');
        textstr = sprintf('%s',event_types{cc}); 
        ht = set_axes_top_text_no_line(ff.hf,gca,textstr,[0.01 -0.01 0 0]); set(ht,'Fontsize',6,'FontWeight','Normal');
        if cc == 1
            ylabel('(\mum)');
%             textstr = sprintf('Excited (Exc)'); ht = set_axes_top_text_no_line(ff.hf,gca,textstr,[-0.01 -0.05 0 0]); set(ht,'Fontsize',6,'FontWeight','Normal');
%             textstr = sprintf('Response to Air Onset (Brake Configuration'); ht = set_axes_top_text_no_line(ff.hf,gca,textstr,[-0.01 0.01 0 0]); set(ht,'Fontsize',6,'FontWeight','Bold');
        else
%           textstr = sprintf('Inhibited (Inh)'); ht = set_axes_top_text_no_line(ff.hf,gca,textstr,[-0.01 -0.05 0 0]); set(ht,'Fontsize',6,'FontWeight','Normal');
        end
    end

    save_pdf(ff.hf,mData.pdf_folder,'rasters.pdf',600);

%% For 
an = 1;
si = Lb ;
Rs = o.Rs(:,si);
props1 = get_props_Rs(Rs,ntrials);
R = Rs{an};
sel_pop = cell_list_op(props1,{'inh'});
plotRasters_simplest(R,find(sel_pop{an}))

%% Time rasters

    magfac = mData.magfac;
    ff = makeFigureRowsCols(108,[10 3 1.9 1.3],'RowsCols',[1 2],'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.02 0.20],'widthHeightAdjustment',[10 -580]);
    MY = 10; ysp = 0.5; mY = 0; titletxt = ''; ylabeltxt = {'Trial #'};
    stp = 0.27*magfac; widths = ([1.3 1.3 1.3 1.3 1.3 0.5 0.5 0.5]-0.71)*magfac; gap = 0.24*magfac;
    adjust_axes(ff,[mY MY],stp,widths,gap,{''});

    cellN = [168 352 45]; %21 51
%     an = 2; cn = 3;R = Rs_A{an,cn};cellN = [15 191 72 155 95 64 ];
    cbar_p_shift = [0.01 0.09 -0.05 -0.3];
    for cc = 1:2
        c = cellN(cc);
        thisRaster = R.sp_rasters(:,:,c);

        ax = ff.h_axes(1,cc);
        axes(ax); xlabel(ax,'Time (s)');
        m = min(thisRaster(:));
        M = max(thisRaster(:));
        xdata = [-2 0 2]; ydata = [1 10];
        imagesc(xdata,ydata,thisRaster,[m M]);
        hold on;
        set(gca,'Ydir','normal');
        if cc > 1
                set(gca,'ytick',[]);
        else
            ylabel('Trial #');
        end
        xlabel('Time (s)'); 
        format_axes(gca)
        colormap parula;
        box off;
        ylims = ylim;
        plot([0 0],ylims,'m');
        %****** color bar
        mM = round(min(thisRaster(:)),1); MM = round(max(thisRaster(:)),1);
        [hc,hca] = putColorBar(gca,cbar_p_shift,[mM MM],5,'eastoutside',[0.1 0.11 0.09 0.16]);

        %******* make new axes and plot mean and gaussian fitting
        pos = get(ax,'Position'); ha = axes; set(ha,'units','inches');set(ha,'Position',pos);
        changePosition(ha,[0 pos(4)+0.05 0 -0.35]);

        xs = linspace(-2,2,size(thisRaster,2));
        mSig = nanmean(thisRaster);
        plot(xs,mSig,'b');hold on;
        fitplot = gauss_fit(1:length(xs),R.gauss_fit_on_mean.coefficients_Rs_mean(c,1:3),R.gauss_fit_on_mean.gauss1Formula);
%         plot(xs,fitplot,'linewidth',0.5,'color','m');
        box off;
        ylim([0 round(max(mSig)+(max(mSig)/10),1)])
        set(ha,'xtick',[],'ytick',round(max(mSig),1));
        format_axes(gca);
        ylims = ylim;
        plot([0 0],ylims,'m');
        textstr = sprintf('Cell %d (%.1f, %.1f)',c,R.info_metrics.ShannonMI_Zsh(c),R.gauss_fit_on_mean.coefficients_Rs_mean(c,4));
        textstr = sprintf('Cell %d (%.1f)',c,R.info_metrics.ShannonMI_Zsh(c));
%         textstr = sprintf('Cell %d',c); ht = set_axes_top_text_no_line(ff.hf,gca,textstr,[-0.01 -0.05 0 0]); set(ht,'Fontsize',6,'FontWeight','Normal');
        if cc == 1
            ylabel('FR (AU)');
            textstr = sprintf('Excited (Exc)'); ht = set_axes_top_text_no_line(ff.hf,gca,textstr,[-0.01 -0.05 0 0]); set(ht,'Fontsize',6,'FontWeight','Normal');
%             textstr = sprintf('Response to Air Onset (Brake Configuration'); ht = set_axes_top_text_no_line(ff.hf,gca,textstr,[-0.01 0.01 0 0]); set(ht,'Fontsize',6,'FontWeight','Bold');
        else
          textstr = sprintf('Inhibited (Inh)'); ht = set_axes_top_text_no_line(ff.hf,gca,textstr,[-0.01 -0.05 0 0]); set(ht,'Fontsize',6,'FontWeight','Normal');
        end
    end

    save_pdf(ff.hf,mData.pdf_folder,'rasters.pdf',600);
    
    
%% For 
an = 4;
si = M_On ;
Rs = o.Rs(:,si);
props1 = get_props_Rs(Rs,50);
R = Rs{an};
sel_pop = cell_list_op(props1,{'exc'});
plotRasters_simplest(R,find(sel_pop{an}))

%% Time rasters

    magfac = mData.magfac;
    ff = makeFigureRowsCols(108,[10 3 1.9 1.3],'RowsCols',[1 2],'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.02 0.20],'widthHeightAdjustment',[10 -580]);
    MY = 10; ysp = 0.5; mY = 0; titletxt = ''; ylabeltxt = {'Trial #'};
    stp = 0.27*magfac; widths = ([1.3 1.3 1.3 1.3 1.3 0.5 0.5 0.5]-0.71)*magfac; gap = 0.24*magfac;
    adjust_axes(ff,[mY MY],stp,widths,gap,{''});

    cellN = [97 57]; %21 51
%     an = 2; cn = 3;R = Rs_A{an,cn};cellN = [15 191 72 155 95 64 ];
    cbar_p_shift = [0.01 0.09 -0.05 -0.3];
    for cc = 1:2
        c = cellN(cc);
        thisRaster = R.sp_rasters(:,:,c);

        ax = ff.h_axes(1,cc);
        axes(ax); xlabel(ax,'Time (s)');
        m = min(thisRaster(:));
        M = max(thisRaster(:));
        xdata = [-1.5 0 1.5]; ydata = [1 10];
        imagesc(xdata,ydata,thisRaster,[m M]);
        hold on;
        set(gca,'Ydir','normal');
        if cc > 1
                set(gca,'ytick',[]);
        else
            ylabel('Trial #');
        end
        xlabel('Time (s)'); 
        format_axes(gca)
        colormap parula;
        box off;
        ylims = ylim;
        set(gca,'xtick',[-1.5 0 1.5],'xticklabel',{'-1.5','0','1.5'});
        plot([0 0],ylims,'m');
        %****** color bar
        mM = round(min(thisRaster(:)),1); MM = round(max(thisRaster(:)),1);
        [hc,hca] = putColorBar(gca,cbar_p_shift,[mM MM],5,'eastoutside',[0.1 0.11 0.09 0.16]);

        %******* make new axes and plot mean and gaussian fitting
        pos = get(ax,'Position'); ha = axes; set(ha,'units','inches');set(ha,'Position',pos);
        changePosition(ha,[0 pos(4)+0.05 0 -0.35]);

        xs = linspace(-1.5,1.5,size(thisRaster,2));
        mSig = nanmean(thisRaster);
        plot(xs,mSig,'b');hold on;
        fitplot = gauss_fit(1:length(xs),R.gauss_fit_on_mean.coefficients_Rs_mean(c,1:3),R.gauss_fit_on_mean.gauss1Formula);
%         plot(xs,fitplot,'linewidth',0.5,'color','m');
        box off;
        ylim([0 round(max(mSig)+(max(mSig)/10),1)])
        set(ha,'xtick',[],'ytick',round(max(mSig),1));
        format_axes(gca);
        ylims = ylim;
        plot([0 0],ylims,'m');
        textstr = sprintf('Cell %d (%.1f, %.1f)',c,R.info_metrics.ShannonMI_Zsh(c),R.gauss_fit_on_mean.coefficients_Rs_mean(c,4));
        textstr = sprintf('Cell %d (%.1f)',c,R.info_metrics.ShannonMI_Zsh(c));
%         textstr = sprintf('Cell %d',c); ht = set_axes_top_text_no_line(ff.hf,gca,textstr,[-0.01 -0.05 0 0]); set(ht,'Fontsize',6,'FontWeight','Normal');
        if cc == 1
            ylabel('FR (AU)');
            textstr = sprintf('Excited (Exc)'); ht = set_axes_top_text_no_line(ff.hf,gca,textstr,[-0.01 -0.05 0 0]); set(ht,'Fontsize',6,'FontWeight','Normal');
%             textstr = sprintf('Response to Air Onset (Brake Configuration'); ht = set_axes_top_text_no_line(ff.hf,gca,textstr,[-0.01 0.01 0 0]); set(ht,'Fontsize',6,'FontWeight','Bold');
        else
          textstr = sprintf('Inhibited (Inh)'); ht = set_axes_top_text_no_line(ff.hf,gca,textstr,[-0.01 -0.05 0 0]); set(ht,'Fontsize',6,'FontWeight','Normal');
        end
        xlim([min(xs) max(xs)]);
    end

    save_pdf(ff.hf,mData.pdf_folder,'rasters.pdf',600);
