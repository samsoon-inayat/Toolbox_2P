function raster_properties

%% general for all properties including responsivity, response fidelity, zMI, Rs
ntrials = 50;
si = [C1_t_D C2_t_D C3_t_D C4_t_D];
% si = [C1_i_T C2_i_T C3_i_T C4_i_T];
Rs_C = oC.Rs(:,si); Rs_A = oA.Rs(:,si); mRs_C = oC.mR(:,si); mRs_A = oA.mR(:,si);
props_C = get_props_Rs(oC.Rs(:,si),ntrials); props_A = get_props_Rs(oA.Rs(:,si),ntrials);
% pop_var_name = {'all','vals','valsT','Nvals','good_zMI','Ngood_zMI'};
pop_var_name = {'good_zMI','good_Gauss','good_MFR'};
pop_var_name = {'all'};
pop_var_name = {'vals'};
% pop_var_name = {'vals','good_zMI'};
% pop_var_name = {'vals','good_zMI','good_Gauss'};
sel_pop_C = cell_list_op(props_C,pop_var_name); sel_pop_A = cell_list_op(props_A,pop_var_name);
% pop_var_name = {'Nvals'};
% sel_pop_CNR = cell_list_op(props_C,pop_var_name); sel_pop_ANR = cell_list_op(props_A,pop_var_name);

params = {'perc','N_Resp_Trials','zMI','rs','nan_zMI','nan_rs','HaFD','HiFD','PWs','centers','peak_locations','mean_FR','MFR'};
params = {'perc','N_Resp_Trials','zMI','rs','PWs','centers','peak_locations','mean_FR','MFR'};
varT = 3;%:length(params)
[~,~,pop_C] = plotDistributions(sel_pop_C);  [~,~,pop_A] = plotDistributions(sel_pop_A);
eval(sprintf('var_CT = props_C.%s;',params{varT}));  eval(sprintf('var_AT = props_A.%s;',params{varT}));
% [~,~,var_C] = plotDistributions(var_CT);  [~,~,var_A] = plotDistributions(var_AT);
out_C = get_vals(var_CT,sel_pop_C); out_A = get_vals(var_AT,sel_pop_A); 

xlabels = {NaN,{'RF (%)'},'zMI','R-Sq','FW (cm)','Cen (cm)','PL (cm)'};
ylabels = {NaN,{'RF (%)'},'zMI','R-Sq','FW (cm)','Cen (cm)','PL (cm)'};
incrs = [NaN,10,0.51,0.01,1,15,1];
mYs = [NaN,0,-0.5,0,0,0,0];
MYs = [NaN,90,0.8125,0.75,2,20,20];
ysps = [NaN,10,0.15,0.2,0.15,0.15,3];

magfac = mData.magfac;
ff = makeFigureRowsCols(108,[4 5 3.45 1.0],'RowsCols',[1 4],'spaceRowsCols',[0.03 -0.01],'rightUpShifts',[0.1 0.25],'widthHeightAdjustment',[10 -500]);
MY = 8; ysp = 1; mY = 0; 
stp = 0.35*magfac; widths = ([0.5 0.5 0.5 0.5 0.25 0.25 0.25 0.25]+0.12)*magfac; gap = 0.175*magfac;
adjust_axes(ff,[mY MY],stp,widths,gap,{'Cell #'});
% shift_axes(ff,5:8,0.35,gap);
var_CvD = out_C(:,1:4); var_AvD = out_A(:,1:4);
for ci = 1:4
    distD = [var_CvD(:,ci) var_AvD(:,ci)];
%     [~,~,var_Ctt] = plotDistributions(var_CvD(:,ci)); [~,~,var_Att] = plotDistributions(var_AvD(:,ci));
%     [h,p,ks2stat] = kstest2(var_Ctt{1},var_Att{1});
    tcolors = {'k','r'};
    [distDo,allVals,allValsG] = plotDistributions(distD);
    minBin = min(allVals);
    maxBin = max(allVals);
    incr = incrs(varT);
%         hf = figure(8);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.85 1],'color','w');
    axes(ff.h_axes(1,ci)); hold on;
    [ha,hb,hca,ylims] = plotDistributions(allValsG,'colors',tcolors,'maxY',0.5,'min',minBin,'incr',incr,'max',maxBin,'do_mean','No');
%     [ha,hb,hca] = plotDistributions(distDo,'colors',tcolors,'maxY',100,'min',minBin,'incr',incr,'max',maxBin,'do_mean','Yes');
    set(gca,'FontSize',6,'FontWeight','Bold','TickDir','out','xcolor','k','ycolor','k');
%     changePosition(gca,[0.129 0.15 -0.09 -0.13]);
    ylims
    ylim([0 100]); xlim([minBin maxBin]);
    if ci == 1
        put_axes_labels(ha,{xlabels{varT},[0 0 0]},{{'Cells (%)'},[0 0 0]});
    else
%         put_axes_labels(ha,{'Firing Rate (A.U.)',[0 0 0]},{{'Neurons (%)'},[0 0 0]});
        put_axes_labels(ha,{xlabels{varT},[0 0 0]},{{''},[0 0 0]});
    end
    format_axes_b(ha);
    [ks2.h,ks2.p,ks2.ks2stat] = kstest2(allValsG{1},allValsG{2}); ks2.DF1 = length(allValsG{1}); ks2.DF2 = length(allValsG{2});
    print_for_manuscript(ks2,'KS2');
    [t2.h,t2.p,t2coi,t2.tstat] = ttest2(allValsG{1},allValsG{2}); t2.cd = computeCohen_d(allValsG{1},allValsG{2});
    print_for_manuscript(t2,'t2');
    
    mVar(1,1) = mean(allValsG{1}); mVar(2,1) = mean(allValsG{2});
    semVar(1,1) = std(allValsG{1})/sqrt(length(allValsG{1})); semVar(2,1) = std(allValsG{2})/sqrt(length(allValsG{2}));
    ylims = ylim;
%     plot([mVar(1) mVar(1)],ylims,'color',tcolors{1});
%     plot([mVar(2) mVar(2)],ylims,'color',tcolors{2});
%     ht = set_axes_top_text_no_line(gcf,ha,'KS-Test',[0.0 -0.01 0 0]);set(ht,'FontSize',7);
    titletxt = sprintf('K, %s (T, %s)',getNumberOfAsterisks(ks2.p),getNumberOfAsterisks(t2.p));
    ht = set_axes_top_text_no_line(gcf,ha,titletxt,[0.0 -0.01 0 0]);set(ht,'FontSize',7);
    titletxt = sprintf('C%d',ci);
    ht = set_axes_top_text_no_line(gcf,ha,titletxt,[0.031 0.1 0 0]);set(ht,'FontSize',8);
    if 1%varT == 2 
        titletxt = sprintf('n = %d,',length(allValsG{1}));
        ht = set_axes_top_text_no_line(gcf,ha,titletxt,[0.03 -0.5 0 0]);set(ht,'FontSize',6,'Color','k');
        titletxt = sprintf('%d',length(allValsG{2}));
        ht = set_axes_top_text_no_line(gcf,ha,titletxt,[0.125 -0.5 0 0]);set(ht,'FontSize',6,'Color','r');
    end
%     if varT == 5 
%         titletxt = sprintf('n = %d,',length(allValsG{1}));
%         ht = set_axes_top_text_no_line(gcf,ha,titletxt,[0.005 -0.41 0 0]);set(ht,'FontSize',7,'Color','k');
%         titletxt = sprintf('%d',length(allValsG{2}));
%         ht = set_axes_top_text_no_line(gcf,ha,titletxt,[0.055 -0.41 0 0]);set(ht,'FontSize',7,'Color','r');
%     end
end
save_pdf(ff.hf,mData.pdf_folder,sprintf('Distribution_firing_rate'),600);
