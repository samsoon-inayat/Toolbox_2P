function firing_rate_motion_vs_rest

% mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size;
% ei_C = evalin('base','ei10_C'); 
% ei_A = evalin('base','ei10_A'); 

% selContexts = [1 2 3 4];
% rasterNames = {'airD','airD','airD','airD'};

Rs_C = oC.Rs;% get_rasters_data(ei_C,selContexts,rasterNames);
Rs_A = oA.Rs;% get_rasters_data(ei_A,selContexts,rasterNames);
% typeP = {'all','vals'
typeP  = 'all';
thr = -1;
%%
ntrials = 50;
si = [C1_t_D C2_t_D C3_t_D C4_t_D C1_i_T C2_i_T C3_i_T C4_i_T];
Rs_C = oC.Rs(:,si); Rs_A = oA.Rs(:,si); mRs_C = oC.mR(:,si); mRs_A = oA.mR(:,si);
props_C = get_props_Rs(oC.Rs(:,si),ntrials); props_A = get_props_Rs(oA.Rs(:,si),ntrials);
% pop_var_name = {'all','vals','valsT','Nvals','good_zMI','Ngood_zMI'};
pop_var_name = {'good_zMI','good_Gauss','good_MFR'};
% pop_var_name = {'all'};
pop_var_name = {'vals'};
sel_pop_C = cell_list_op(props_C,pop_var_name); sel_pop_A = cell_list_op(props_A,pop_var_name);
%%
out_C = get_spike_rate_ext(ei_C,thr,sel_pop_C);
out_A = get_spike_rate_ext(ei_A,thr,sel_pop_A);
tcolors = {'k','r','k','r'};
n=0;
%% 
magfac = mData.magfac;
ff = makeFigureRowsCols(108,[6 3 6.9 1.25],'RowsCols',[1 4],'spaceRowsCols',[0.03 -0.01],'rightUpShifts',[0.1 0.275],'widthHeightAdjustment',[10 -385]);
MY = 8; ysp = 1; mY = 0; 
stp = 0.45*magfac; widths = ([0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5]+0.33)*magfac; gap = 0.79*magfac;
adjust_axes(ff,[mY MY],stp,widths,gap,{'Cell #'});
var_CvD = out_C; var_AvD = out_A; 
for ci = 1:4
    distD = [var_CvD(:,ci) var_AvD(:,ci)];
    [~,~,var_Ctt] = plotDistributions(var_CvD(:,ci)); [~,~,var_Att] = plotDistributions(var_AvD(:,ci));
    [h,p,ks2stat] = kstest2(var_Ctt{1},var_Att{1});
    tcolors = {'k','r'};
    [distDo,allVals,allValsG] = plotDistributions(distD);
    minBin = min(allVals);
    maxBin = max(allVals);
    incr = 0.001;
%         hf = figure(8);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.85 1],'color','w');
    axes(ff.h_axes(1,ci)); hold on;
%    [ha,hb,hca] = plotDistributions(distD,'colors',tcolors,'maxY',maxBin,'cumPos',[0.5 0.26 0.25 0.5],'min',minBin,'incr',incr,'max',maxBin);
    [ha,hb,hca] = plotDistributions(distDo,'colors',tcolors,'maxY',100,'min',minBin,'incr',incr,'max',maxBin);
    set(gca,'FontSize',6,'FontWeight','Bold','TickDir','out','xcolor','k','ycolor','k');
%     changePosition(gca,[0.129 0.15 -0.09 -0.13]);
    ylim([0 100]);
    if ci == 1
        put_axes_labels(ha,{'Firing Rate (A.U.)',[0 0 0]},{{'Neurons (%)'},[0 0 0]});
    else
        put_axes_labels(ha,{'Firing Rate (A.U.)',[0 0 0]},{{''},[0 0 0]});
    end
    format_axes_b(ha);
    xlim([minBin maxBin]);
    pos = get(ha,'Position'); han = axes; set(han,'units','inches');set(han,'Position',pos);
    changePosition(han,[widths(1)+0.3 0 -0.5 0]);
    box off;
    format_axes_b(han);
    [h,p,ks2stat] = kstest2(allValsG{1},allValsG{2});
    titletxt = sprintf('%.3f',p);
%     titletxt = [titletxt 
    ht = set_axes_top_text_no_line(gcf,ha,titletxt,[0 -0.051 0 0]);
end
save_pdf(ff.hf,mData.pdf_folder,sprintf('Distribution_firing_rate'),600);