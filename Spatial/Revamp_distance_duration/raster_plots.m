function raster_plots

%% rasters speed
magfac = mData.magfac;
ff = makeFigureRowsCols(108,[10 3 6.99 2],'RowsCols',[2 6],'spaceRowsCols',[0.3 0.02],'rightUpShifts',[0.02 0.135],'widthHeightAdjustment',[10 -325]);
MY = 10; ysp = 0.5; mY = 0; titletxt = ''; ylabeltxt = {'Trial #'};
stp = 0.35*magfac; widths = (ones(1,18)*0.77)*magfac; gap = 0.32*magfac;
adjust_axes(ff,[mY MY],stp,widths,gap,{''});
shift_axes(ff,[4 5 6;4 5 6],0.2,gap);
pos = get(ff.h_axes(1,1),'Position');
pos_msig = [0 pos(4)+0.05 0 -0.15];
cbar_t_shift = [0.1 0.2 0.09 0.2];
%  time
ntrials = 50; 
si = [Ar_t_T ArL_t_T Ars_t_T Ar_i_T ArL_i_T Ars_i_T];
Rs_C = o.Rs(:,si);mRs_C = o.mR(:,si);
props_C = get_props_Rs(Rs_C,ntrials);
pop_var_name = {'all','vals','valsT','Nvals','good_zMI','Ngood_zMI'};
pop_var_name = {'vals','good_zMI'};
sel_pop_C = cell_list_op(props_C,pop_var_name); 

an = 4; Rs = Rs_C(an,:);
sel_pop = sel_pop_C(an,5); sel_pop_or = cell_list_op(sel_pop,[],'or',1); sel_pop_or = sel_pop_or{1};
cellN = find(sel_pop_or); ci = 12;
cbar_p_shift = [0.01 0.09 -0.05 -0.2];
configN = [3 4 5 3 4 5];
for cc = 1:6
    c = cellN(ci);
    R = Rs{cc};
    thisRaster = R.speed;%sp_rasters(:,:,c);
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
        set(gca,'ytick',[1 5 10]);
        ylabel('Trial #');
    end
    xlabel('Time (s)'); 
    format_axes(gca)
    colormap parula;
    box off;

    %****** color bar
    mM = round(min(thisRaster(:)),1); MM = round(max(thisRaster(:)),1);
    try
    [hc,hca] = putColorBar(gca,cbar_p_shift,[mM MM],5,'eastoutside',cbar_t_shift);
    catch
      continue;
    end

    %******* make new axes and plot mean and gaussian fitting
    pos = get(ax,'Position'); ha = axes; set(ha,'units','inches');set(ha,'Position',pos);
    changePosition(ha,pos_msig);

    xs = linspace(0,sz*bw,sz);
    mSig = nanmean(thisRaster);
    plot(xs,mSig,'b');hold on;

    fitplot = gauss_fit(1:length(xs),R.gauss_fit_on_mean.coefficients_Rs_mean(c,1:3),R.gauss_fit_on_mean.gauss1Formula);
%     plot(xs,fitplot,'linewidth',0.5,'color','m');
    box off;
    try
    ylim([0 round(max(mSig)+(max(mSig)/10),1)])
    catch
      continue;
    end
    set(ha,'xtick',[],'ytick',round(max(mSig),1));
    format_axes(gca);
    textstr = sprintf('Cell %d (%.1f, %.1f)',c,R.info_metrics.ShannonMI_Zsh(c),R.gauss_fit_on_mean.coefficients_Rs_mean(c,4));
    if cc < 4
      textstr = sprintf('C%d - AOn',configN(cc));
    else
      textstr = sprintf('C%d - AOff',configN(cc));
    end
    ht = set_axes_top_text_no_line(ff.hf,gca,textstr,[0.01 -0.09 0 0]); set(ht,'Fontsize',6,'FontWeight','Bold');
    if cc == 1
        ylabel('cm/s');
    end
end

% dist
ntrials = 50; 
si = [Ar_t_D ArL_t_D Ars_t_D Ar_i_D ArL_i_D Ars_i_D];
Rs_C = o.Rs(:,si);mRs_C = o.mR(:,si);
props_C = get_props_Rs(Rs_C,ntrials);
pop_var_name = {'all','vals','valsT','Nvals','good_zMI','Ngood_zMI'};
pop_var_name = {'vals','good_zMI'};
sel_pop_C = cell_list_op(props_C,pop_var_name); 

Rs = Rs_C(an,:);
% sel_pop = sel_pop_C(an,:); sel_pop_or = cell_list_op(sel_pop,[],'or',1); sel_pop_or = sel_pop_or{1};
% cellN = find(sel_pop_or);
% cbar_p_shift = [0.01 0.09 -0.05 -0.3];
for cc = 1:6
    c = cellN(ci);
    R = Rs{cc};
    thisRaster = R.speed;%sp_rasters(:,:,c);
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
        set(gca,'ytick',[1 5 10]);
        ylabel('Trial #');
    end
    xlabel('Distance (cm)'); 
    format_axes(gca)
    colormap parula;
    box off;

    %****** color bar
    mM = round(min(thisRaster(:)),1); MM = round(max(thisRaster(:)),1);
    try
    [hc,hca] = putColorBar(gca,cbar_p_shift,[mM MM],5,'eastoutside',cbar_t_shift);
    catch
      continue;
    end

    %******* make new axes and plot mean and gaussian fitting
    pos = get(ax,'Position'); ha = axes; set(ha,'units','inches');set(ha,'Position',pos);
    changePosition(ha,pos_msig);

    xs = linspace(0,sz*bw,sz);
    mSig = nanmean(thisRaster);
    plot(xs,mSig,'b');hold on;
    fitplot = gauss_fit(1:length(xs),R.gauss_fit_on_mean.coefficients_Rs_mean(c,1:3),R.gauss_fit_on_mean.gauss1Formula);
%     plot(xs,fitplot,'linewidth',0.5,'color','m');
    box off;
    try
    ylim([0 round(max(mSig)+(max(mSig)/10),1)])
    catch
      continue;
    end
    set(ha,'xtick',[],'ytick',round(max(mSig),1));
    format_axes(gca);
%     textstr = sprintf('Cell %d (%.1f, %.1f)',c,R.info_metrics.ShannonMI_Zsh(c),R.gauss_fit_on_mean.coefficients_Rs_mean(c,4));
    textstr = sprintf('%.1f',R.info_metrics.ShannonMI_Zsh(c));
%      ht = set_axes_top_text_no_line(ff.hf,gca,textstr,[0.01 -0.05 0 0]); set(ht,'Fontsize',6,'FontWeight','Bold');
    if cc == 1
        ylabel('cm/s');
    end
end

% save_pdf(ff.hf,mData.pdf_folder,'rasters.pdf',600);



save_pdf(ff.hf,mData.pdf_folder,'rasters.pdf',600);


%% rasters
magfac = mData.magfac;
ff = makeFigureRowsCols(108,[10 3 6.98 2],'RowsCols',[2 6],'spaceRowsCols',[0.3 0.02],'rightUpShifts',[0.02 0.15],'widthHeightAdjustment',[10 -325]);
MY = 10; ysp = 0.5; mY = 0; titletxt = ''; ylabeltxt = {'Trial #'};
stp = 0.35*magfac; widths = (ones(1,18)*0.77)*magfac; gap = 0.32*magfac;
adjust_axes(ff,[mY MY],stp,widths,gap,{''});
shift_axes(ff,[4 5 6;4 5 6],0.2,gap);
pos = get(ff.h_axes(1,1),'Position');
pos_msig = [0 pos(4)+0.05 0 -0.15];
cbar_t_shift = [0.1 0.2 0.09 0.2];
%  time
ntrials = 50; 
si = [Ar_t_T ArL_t_T Ars_t_T Ar_i_T ArL_i_T Ars_i_T];
Rs_C = o.Rs(:,si);mRs_C = o.mR(:,si);
props_C = get_props_Rs(Rs_C,ntrials);
pop_var_name = {'all','vals','valsT','Nvals','good_zMI','Ngood_zMI'};
pop_var_name = {'vals','good_zMI'};
sel_pop_C = cell_list_op(props_C,pop_var_name); 

si_cn_ap_H = si_cn_ap(:,[1 3 5 2 4 6]);

an = 4; Rs = Rs_C(an,:);
sel_pop = sel_pop_C(an,5); sel_pop_or = cell_list_op(sel_pop,[],'or',1); sel_pop_or = sel_pop_or{1};
sel_popD = sel_pop_CD(an,5); sel_pop_orD = cell_list_op(sel_popD,[],'or',1); sel_pop_orD = sel_pop_orD{1};
cellN = find(sel_pop_orD); ci = 28;%19; % 12 17;
cbar_p_shift = [0.01 0.09 -0.05 -0.2];
for cc = 1:6
    c = cellN(ci);
    R = Rs{cc}; MVT = met_valsT{1}{an,si_cn_ap_H(1,cc),si_cn_ap_H(2,cc)};
    zval_here = [MVT.PC(c,2) MVT.MI(c,2)];
    thisRaster = R.sp_rasters(:,:,c);
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
      set(gca,'ytick',[1 5 10]);
        ylabel('Trial #');
    end
    xlabel('Time (s)'); 
    format_axes(gca)
    colormap parula;
    box off;

    %****** color bar
    mM = round(min(thisRaster(:)),1); MM = round(max(thisRaster(:)),1);
    try
    [hc,hca] = putColorBar(gca,cbar_p_shift,[mM MM],5,'eastoutside',cbar_t_shift);
    catch
      continue;
    end

    %******* make new axes and plot mean and gaussian fitting
    pos = get(ax,'Position'); ha = axes; set(ha,'units','inches');set(ha,'Position',pos);
    changePosition(ha,pos_msig);

    xs = linspace(0,sz*bw,sz);
    mSig = nanmean(thisRaster);
    plot(xs,mSig,'b');hold on;

    fitplot = gauss_fit(1:length(xs),R.gauss_fit_on_mean.coefficients_Rs_mean(c,1:3),R.gauss_fit_on_mean.gauss1Formula);
    % plot(xs,fitplot,'linewidth',0.5,'color','m');
    box off;
    YM = round(max(mSig)+(max(mSig)/10),1);
    if YM == 0
      YM = round(max(mSig)+(max(mSig)/10),2);
    end
    ylim([0 YM]);
%     set(ha,'xtick',[],'ytick',round(max(mSig),1));
    set(ha,'xtick',[],'ytick',YM);
    format_axes(gca);
    textstr = sprintf('Cell %d (%.1f, %.1f)',c,R.info_metrics.ShannonMI_Zsh(c),R.gauss_fit_on_mean.coefficients_Rs_mean(c,4));
    % if cc == 1
    %   textstr = sprintf(' %.1f   A%d-C%d',R.info_metrics.ShannonMI_Zsh(c),an,c);
    %   ht = set_axes_top_text_no_line(ff.hf,gca,textstr,[0.0 0.04 0 0]); set(ht,'Fontsize',6,'FontWeight','Bold');
    % else
      textstr = sprintf('%.1f, %.1f',zval_here);
      ht = set_axes_top_text_no_line(ff.hf,gca,textstr,[0.01 0.04 0 0]); set(ht,'Fontsize',6,'FontWeight','Bold');
    % end
    if cc == 1
        ylabel('FR (AU)');
    end
end

% dist
ntrials = 50; 
si = [Ar_t_D ArL_t_D Ars_t_D Ar_i_D ArL_i_D Ars_i_D];
Rs_C = o.Rs(:,si);mRs_C = o.mR(:,si);
% props_C = get_props_Rs(Rs_C,ntrials);
% pop_var_name = {'all','vals','valsT','Nvals','good_zMI','Ngood_zMI'};
% pop_var_name = {'vals','good_zMI'};
% sel_pop_C = cell_list_op(props_C,pop_var_name); 

Rs = Rs_C(an,:);
% sel_pop = sel_pop_C(an,:); sel_pop_or = cell_list_op(sel_pop,[],'or',1); sel_pop_or = sel_pop_or{1};
% cellN = find(sel_pop_or);
% cbar_p_shift = [0.01 0.09 -0.05 -0.3];
for cc = 1:6
    c = cellN(ci);
    R = Rs{cc};MVD = met_valsD{2}{an,si_cn_ap_H(1,cc),si_cn_ap_H(2,cc)};
    zval_here = [MVD.PC(c,2) MVD.MI(c,2)];
    thisRaster = R.sp_rasters(:,:,c);
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
      set(gca,'ytick',[1 5 10]);
        ylabel('Trial #');
    end
    xlabel('Distance (cm)'); 
    format_axes(gca)
    colormap parula;
    box off;

    %****** color bar
    mM = round(min(thisRaster(:)),1); MM = round(max(thisRaster(:)),1);
    try
    [hc,hca] = putColorBar(gca,cbar_p_shift,[mM MM],5,'eastoutside',cbar_t_shift);
    catch
      continue;
    end

    %******* make new axes and plot mean and gaussian fitting
    pos = get(ax,'Position'); ha = axes; set(ha,'units','inches');set(ha,'Position',pos);
    changePosition(ha,pos_msig);

    xs = linspace(0,sz*bw,sz);
    mSig = nanmean(thisRaster);
    plot(xs,mSig,'b');hold on;
    fitplot = gauss_fit(1:length(xs),R.gauss_fit_on_mean.coefficients_Rs_mean(c,1:3),R.gauss_fit_on_mean.gauss1Formula);
    % plot(xs,fitplot,'linewidth',0.5,'color','m');
    box off;
    YM = round(max(mSig)+(max(mSig)/10),1);
    if YM == 0
      YM = 0.1;
    end
    ylim([0 YM])
    set(ha,'xtick',[],'ytick',round(max(mSig),1));
    format_axes(gca);
%     textstr = sprintf('Cell %d (%.1f, %.1f)',c,R.info_metrics.ShannonMI_Zsh(c),R.gauss_fit_on_mean.coefficients_Rs_mean(c,4));
    textstr = sprintf('%.1f',R.info_metrics.ShannonMI_Zsh(c));
    textstr = sprintf('%.1f, %.1f',zval_here);
    ht = set_axes_top_text_no_line(ff.hf,gca,textstr,[0.01 0.04 0 0]); set(ht,'Fontsize',6,'FontWeight','Bold');
    if cc == 1
        ylabel('FR (AU)');
    end
end

% save_pdf(ff.hf,mData.pdf_folder,'rasters.pdf',600);



save_pdf(ff.hf,mData.pdf_folder,'rasters.pdf',600);

%% to check representative dur and dis cells
magfac = mData.magfac;
ff = makeFigureRowsCols(108,[10 3 6.98 2],'RowsCols',[2 6],'spaceRowsCols',[0.3 0.02],'rightUpShifts',[0.02 0.15],'widthHeightAdjustment',[10 -325]);
MY = 10; ysp = 0.5; mY = 0; titletxt = ''; ylabeltxt = {'Trial #'};
stp = 0.35*magfac; widths = (ones(1,18)*0.77)*magfac; gap = 0.32*magfac;
adjust_axes(ff,[mY MY],stp,widths,gap,{''});
shift_axes(ff,[4 5 6;4 5 6],0.2,gap);
pos = get(ff.h_axes(1,1),'Position');
pos_msig = [0 pos(4)+0.05 0 -0.15];
cbar_t_shift = [0.1 0.2 0.09 0.2];
%  time
ntrials = 50; 
si = [Ar_t_T ArL_t_T Ars_t_T Ar_i_T ArL_i_T Ars_i_T];
Rs_C = o.Rs(:,si);mRs_C = o.mR(:,si);
props_C = get_props_Rs(Rs_C,ntrials);
% pop_var_name = {'all','vals','valsT','Nvals','good_zMI','Ngood_zMI'};
% pop_var_name = {'vals','good_zMI'};
% sel_pop_C = cell_list_op(props_C,pop_var_name); 
sel_pop_C = [dur_cells_T dur_cells_I];


siD = [Ar_t_D ArL_t_D Ars_t_D Ar_i_D ArL_i_D Ars_i_D];
Rs_CD = o.Rs(:,siD);mRs_CD = o.mR(:,siD);
props_CD = get_props_Rs(Rs_CD,ntrials);
% sel_pop_CD = cell_list_op(props_CD,pop_var_name); 
sel_pop_CD = [dis_cells_T dis_cells_I];


an = 4; Rs = Rs_C(an,:);
sel_pop = sel_pop_C(an,4); sel_pop_or = cell_list_op(sel_pop,[],'or',1); sel_pop_or = sel_pop_or{1};
sel_popD = sel_pop_CD(an,1); sel_pop_orD = cell_list_op(sel_popD,[],'or',1); sel_pop_orD = sel_pop_orD{1};
cellN = find(sel_pop_orD); ci = 12;%19; % 12 17;
cbar_p_shift = [0.01 0.09 -0.05 -0.2];
for cc = 1:6
    c = cellN(ci);
    R = Rs{cc};
    thisRaster = R.sp_rasters(:,:,c);
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
      set(gca,'ytick',[1 5 10]);
        ylabel('Trial #');
    end
    xlabel('Time (s)'); 
    format_axes(gca)
    colormap parula;
    box off;

    %****** color bar
    mM = round(min(thisRaster(:)),1); MM = round(max(thisRaster(:)),1);
    try
    [hc,hca] = putColorBar(gca,cbar_p_shift,[mM MM],5,'eastoutside',cbar_t_shift);
    catch
      continue;
    end

    %******* make new axes and plot mean and gaussian fitting
    pos = get(ax,'Position'); ha = axes; set(ha,'units','inches');set(ha,'Position',pos);
    changePosition(ha,pos_msig);

    xs = linspace(0,sz*bw,sz);
    mSig = nanmean(thisRaster);
    plot(xs,mSig,'b');hold on;

    fitplot = gauss_fit(1:length(xs),R.gauss_fit_on_mean.coefficients_Rs_mean(c,1:3),R.gauss_fit_on_mean.gauss1Formula);
    plot(xs,fitplot,'linewidth',0.5,'color','m');
    box off;
    YM = round(max(mSig)+(max(mSig)/10),1);
    if YM == 0
      YM = round(max(mSig)+(max(mSig)/10),2);
    end
    ylim([0 YM]);
%     set(ha,'xtick',[],'ytick',round(max(mSig),1));
    set(ha,'xtick',[],'ytick',YM);
    format_axes(gca);
    textstr = sprintf('Cell %d (%.1f, %.1f)',c,R.info_metrics.ShannonMI_Zsh(c),R.gauss_fit_on_mean.coefficients_Rs_mean(c,4));
    if cc == 1
      textstr = sprintf(' %.1f   A%d-C%d',R.info_metrics.ShannonMI_Zsh(c),an,c);
      ht = set_axes_top_text_no_line(ff.hf,gca,textstr,[0.0 -0.09 0 0]); set(ht,'Fontsize',6,'FontWeight','Bold');
    else
      textstr = sprintf('%.1f',R.info_metrics.ShannonMI_Zsh(c));
      ht = set_axes_top_text_no_line(ff.hf,gca,textstr,[0.01 -0.09 0 0]); set(ht,'Fontsize',6,'FontWeight','Bold');
    end
    if cc == 1
        ylabel('FR (AU)');
    end
end

% dist
ntrials = 50; 
si = [Ar_t_D ArL_t_D Ars_t_D Ar_i_D ArL_i_D Ars_i_D];
Rs_C = o.Rs(:,si);mRs_C = o.mR(:,si);
% props_C = get_props_Rs(Rs_C,ntrials);
% pop_var_name = {'all','vals','valsT','Nvals','good_zMI','Ngood_zMI'};
% pop_var_name = {'vals','good_zMI'};
% sel_pop_C = cell_list_op(props_C,pop_var_name); 

Rs = Rs_C(an,:);
% sel_pop = sel_pop_C(an,:); sel_pop_or = cell_list_op(sel_pop,[],'or',1); sel_pop_or = sel_pop_or{1};
% cellN = find(sel_pop_or);
% cbar_p_shift = [0.01 0.09 -0.05 -0.3];
for cc = 1:6
    c = cellN(ci);
    R = Rs{cc};
    thisRaster = R.sp_rasters(:,:,c);
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
      set(gca,'ytick',[1 5 10]);
        ylabel('Trial #');
    end
    xlabel('Distance (cm)'); 
    format_axes(gca)
    colormap parula;
    box off;

    %****** color bar
    mM = round(min(thisRaster(:)),1); MM = round(max(thisRaster(:)),1);
    try
    [hc,hca] = putColorBar(gca,cbar_p_shift,[mM MM],5,'eastoutside',cbar_t_shift);
    catch
      continue;
    end

    %******* make new axes and plot mean and gaussian fitting
    pos = get(ax,'Position'); ha = axes; set(ha,'units','inches');set(ha,'Position',pos);
    changePosition(ha,pos_msig);

    xs = linspace(0,sz*bw,sz);
    mSig = nanmean(thisRaster);
    plot(xs,mSig,'b');hold on;
    fitplot = gauss_fit(1:length(xs),R.gauss_fit_on_mean.coefficients_Rs_mean(c,1:3),R.gauss_fit_on_mean.gauss1Formula);
    plot(xs,fitplot,'linewidth',0.5,'color','m');
    box off;
    YM = round(max(mSig)+(max(mSig)/10),1);
    if YM == 0
      YM = 0.1;
    end
    ylim([0 YM])
    set(ha,'xtick',[],'ytick',round(max(mSig),1));
    format_axes(gca);
%     textstr = sprintf('Cell %d (%.1f, %.1f)',c,R.info_metrics.ShannonMI_Zsh(c),R.gauss_fit_on_mean.coefficients_Rs_mean(c,4));
    textstr = sprintf('%.1f',R.info_metrics.ShannonMI_Zsh(c));
    ht = set_axes_top_text_no_line(ff.hf,gca,textstr,[0.01 -0.05 0 0]); set(ht,'Fontsize',6,'FontWeight','Bold');
    if cc == 1
        ylabel('FR (AU)');
    end
end

% save_pdf(ff.hf,mData.pdf_folder,'rasters.pdf',600);



save_pdf(ff.hf,mData.pdf_folder,'rasters.pdf',600);

winopen(fullfile(mData.pdf_folder,'rasters.pdf'));

%% playing with how rasters are made and why is there a difference in response fidelity between distance and time rasters

eei = ei(4);
pp = 1;
onsets = eei{1}.plane{pp}.contexts(5).markers.airI_onsets;
offsets = eei{1}.plane{pp}.contexts(5).markers.airI_offsets;
rasterType = 'time'; binwidths = evalin('base','binwidths');
rastersTS =  make_rasters_quick(eei{1},pp,onsets,offsets,rasterType,binwidths);
rasterType = 'dist'; 
rastersTD =  make_rasters_quick(eei{1},pp,onsets,offsets,rasterType,binwidths);
% conclusion --> there are NaNs in the distance raster for the particular
% case I had been looking (Configuration 5 and Cell 57) ... NaNs are
% because the animal didn't move. It stopped in bin 8 trial 2 that I looked
% at

%% playing with MI of rasters
%% finding dzMIs for trials and intertrials
RsDt = o.Rs(:,[Ar_t_D ArL_t_D Ars_t_D]);  RsTt = o.Rs(:,[Ar_t_T ArL_t_T Ars_t_T]);
RsDi = o.Rs(:,[Ar_i_D ArL_i_D Ars_i_D]);  RsTi = o.Rs(:,[Ar_i_T ArL_i_T Ars_i_T]);
[dzMI_FD,dzMI_FT] = get_zMI_comp_dist_time(RsDt,RsTt,RsDi,RsTi);
%%
an = 4; cn = 1;
tRsD = RsDt{an,cn};

tRsT = RsTt{an,cn};

var1 = tRsT.sp_rasters;
var2 = tRsT.dist;
var3 = tRsT.speed;
var4 = tRsT.space;
%%
for cni = 1:size(var1,3)
    frvar = var1(:,:,cni);
    MIsD(cni) = calculate_MI(frvar,var2,4);
    MIsT(cni) = calculate_MI(frvar,[],4);
    MIsV(cni) = calculate_MI(frvar,var3,4);
    MIsS(cni) = calculate_MI(frvar,var4,4);
end

% MIsD = tRsD.info_metrics.ShannonMI;
%%
figure(1000);clf;hold on;
% plot(MIsT,'r');hold on;
plot(MIsD,'b');
% plot(MIsV,'m');
plot(MIsS,'c');

%% figu
figure(1000);clf;
imagesc(var3);

