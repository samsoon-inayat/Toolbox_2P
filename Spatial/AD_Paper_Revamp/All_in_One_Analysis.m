function All_in_One_Analysis

%% define configurations based on raster plots to analyze and logical vectors for cells
ntrials = 50;
si = [C1_t_D C2_t_D C3_t_D C4_t_D C1_i_T C2_i_T C3_i_T C4_i_T];
Rs_C = oC.Rs(:,si); Rs_A = oA.Rs(:,si); mRs_C = oC.mR(:,si); mRs_A = oA.mR(:,si);
props_C = get_props_Rs(oC.Rs(:,si),ntrials); props_A = get_props_Rs(oA.Rs(:,si),ntrials);
% pop_var_name = {'all','vals','valsT','Nvals','good_zMI','Ngood_zMI'};
pop_var_name = {'good_zMI','good_Gauss','good_MFR'};

% Define cell population
pop_var_name = {'all'}; % all cells
% pop_var_name = {'vals'}; % responsive cells
% pop_var_name = {'vals','good_zMI'}; % highly tuned cells

% get cell ids logical vectors for cell lists
sel_pop_C = cell_list_op(props_C,pop_var_name); sel_pop_A = cell_list_op(props_A,pop_var_name);
%% get properties of selected cells
bins = [0 50 100 150]; bins = [0 150];
nbins = length(bins)-1;
%     sel_pop_C = bin_cells(props_C.peak_locations,bins,sel_pop_C); sel_pop_A = bin_cells(props_A.peak_locations,bins,sel_pop_A);
if nbins > 1
%     [sel_pop_C,perc_C] = bin_cells(props_C.peak_locations,bins,sel_pop_C); [sel_pop_A,perc_A] = bin_cells(props_A.peak_locations,bins,sel_pop_A);
    [sel_pop_C,perc_C] = bin_cells(props_C.centersD,bins,sel_pop_C); [sel_pop_A,perc_A] = bin_cells(props_A.centersD,bins,sel_pop_A);
    ind = 1;
    for ii = 1:4
        for jj = 1:nbins
            temp_tcolors{ind} = mData.colors{ii};
            ind = ind + 1;
        end
    end
end
%     
all_var_C = [];
params = {'perc','N_Resp_Trials','zMI','rs','nan_zMI','nan_rs','HaFD','HiFD','PWs','centers','peak_locations','mean_FR','MFR'};
varT = [2 3 4];%:length(params)
for piii = 1:length(varT)
    pii = varT(piii);
    eval(sprintf('var_C = get_vals(props_C.%s,sel_pop_C);',params{pii})); eval(sprintf('var_A = get_vals(props_A.%s,sel_pop_A);',params{pii}));
    all_var_C{piii} = var_C; all_var_A{piii} = var_A;
end

%% get firing rates and trasients
filename = fullfile(mData.pd_folder,'FR_Transients_All.mat');
% filename = fullfile(mData.pd_folder,'FR_Transients_Responsive.mat');
% filename = fullfile(mData.pd_folder,'FR_Transients_Resp_Highly_tuned.mat');

load(filename); % this should output two variables out_C and out_CT ... same for AD group
%% make scatter and linear regression plots between FR and Transients
% Select Air ON or Air OFF phases
magfac = mData.magfac;
ff = makeFigureRowsCols(108,[4 5 3 1],'RowsCols',[1 8],'spaceRowsCols',[0.03 -0.01],'rightUpShifts',[0.1 0.3],'widthHeightAdjustment',[10 -625]);
MY = 10; ysp = 1; mY = 0; 
stp = 0.4*magfac; widths = ([0.5 0.5 0.5 0.5 0.25 0.25 0.25 0.25]+0.0015)*magfac; gap = 0.175*magfac;
adjust_axes(ff,[mY MY],stp,widths,gap,{''});
shift_axes(ff,5:8,0.35,gap);
confsi = 1:4; %confsi = 5:8;
var_CvD = out_CT(:,confsi); var_AvD = out_AT(:,confsi); var_CvDFR = out_C(:,confsi); var_AvDFR = out_A(:,confsi);
for ci = 1:4
    distD = [var_CvD(:,ci) var_AvD(:,ci)]; distDFR = [var_CvDFR(:,ci) var_AvDFR(:,ci)];
    [distDo,allVals,allValsG] = plotDistributions(distD);  [distDoFR,allValsFR,allValsGFR] = plotDistributions(distDFR);
    trC = allValsG{1}; trA = allValsG{2};  trCFR = allValsGFR{1}; trAFR = allValsGFR{2};
    minBin = min(allVals);
    axes(ff.h_axes(1,ci)); hold on;
    ha = gca;
    sc = scatter(trCFR,trC,1,'.','k');% sc.MarkerEdgeAlpha = 0.5; sc.MarkerFaceAlpha = 0.5
    sc = scatter(trAFR,trA,1,'.','r');% sc.MarkerEdgeAlpha = 0.5; sc.MarkerFaceAlpha = 0.5
    
    opts = statset('Display','iter','TolFun',1e-10);
%     mdl = fitnlm(trCFR,trC,modelfun,beta0,'Options',opts); mdl_C{ci} = mdl;
    grps = [ones(size(trCFR));(2*ones(size(trAFR)))];
    mdlT = table([trCFR;trAFR],[trC;trA],[ones(size(trCFR));(2*ones(size(trAFR)))]);
    mdlT.Properties.VariableNames = {'FR','TR','GR'};
    mdlT.GR = categorical(mdlT.GR);
    fit = fitlm(mdlT,'TR~FR*GR');
    ant = anova(fit);
    pvals(ci) = ant{3,5};
    minFR = min(mdlT{:,1}); maxFR = max(mdlT{:,1});
    xFR = [minFR:0.001:maxFR]'; gFR = ones(size(xFR)); 
    xFR = [xFR;xFR]; gFR = [gFR;(2*gFR)];
    xdT = table(xFR,gFR); xdT.Properties.VariableNames = mdlT.Properties.VariableNames([1 3]); xdT.GR = categorical(xdT.GR);
    trCp = predict(fit,xdT);
    plot(xFR(gFR == 1),trCp(gFR == 1),'k');
    plot(xFR(gFR == 2),trCp(gFR == 2),'r');
%     plot(trAFR,trCp(gFR == 2),'r');
    ylim([0 max(mdlT{:,2})]);
    set(gca,'FontSize',6,'FontWeight','Bold','TickDir','out','xcolor','k','ycolor','k');
    if ci == 1
        put_axes_labels(ha,{{'Avg. FR(A.U.)'},[0 0 0]},{{'Avg. # of Trans./min'},[0 0 0]});
    else
        put_axes_labels(ha,{{'Avg. FR(A.U.)'},[0 0 0]},{{''},[0 0 0]});
    end
    titletxt = sprintf('C%d',ci);
    ht = set_axes_top_text_no_line(gcf,ha,titletxt,[0.031 0.1 0 0]);set(ht,'FontSize',8);
    titletxt = sprintf('%s',getNumberOfAsterisks(pvals(ci))); ht = set_axes_top_text_no_line(gcf,ha,titletxt,[0.031 -0.031 0 0]);set(ht,'FontSize',8);
    format_axes_b(gca);
end
axes(ff.h_axes(1,1));
ht = set_axes_top_text(gcf,ff.h_axes,'Air-ON',[0.051 0.1 0.5 0]);set(ht,'FontSize',8);
save_pdf(ff.hf,mData.pdf_folder,sprintf('Distribution_firing_rate'),600);
%% make scatter and linear regression plots between FR and zMI
% Select Air ON or Air OFF phases
magfac = mData.magfac;
ff = makeFigureRowsCols(108,[4 5 3 1],'RowsCols',[1 8],'spaceRowsCols',[0.03 -0.01],'rightUpShifts',[0.1 0.3],'widthHeightAdjustment',[10 -585]);
MY = 10; ysp = 1; mY = 0; 
stp = 0.4*magfac; widths = ([0.5 0.5 0.5 0.5 0.25 0.25 0.25 0.25]+0.0015)*magfac; gap = 0.175*magfac;
adjust_axes(ff,[mY MY],stp,widths,gap,{''});
shift_axes(ff,5:8,0.35,gap);
confsi = 1:4;% confsi = 5:8
var_CvD = var_C(:,confsi); var_AvD = var_A(:,confsi); var_CvDFR = out_CT(:,confsi); var_AvDFR = out_AT(:,confsi);
for ci = 1:4
    distD = [var_CvD(:,ci) var_AvD(:,ci)]; distDFR = [var_CvDFR(:,ci) var_AvDFR(:,ci)];
    [distDo,allVals,allValsG] = plotDistributions(distD);  [distDoFR,allValsFR,allValsGFR] = plotDistributions(distDFR);
    trC = allValsG{1}; trA = allValsG{2};  trCFR = allValsGFR{1}; trAFR = allValsGFR{2};
    minBin = min(allVals);
    axes(ff.h_axes(1,ci)); hold on;
    ha = gca;
    sc = scatter(trCFR,trC,1,'.','k');% sc.MarkerEdgeAlpha = 0.5; sc.MarkerFaceAlpha = 0.5
    sc = scatter(trAFR,trA,1,'.','r');% sc.MarkerEdgeAlpha = 0.5; sc.MarkerFaceAlpha = 0.5
    
    opts = statset('Display','iter','TolFun',1e-10);
%     mdl = fitnlm(trCFR,trC,modelfun,beta0,'Options',opts); mdl_C{ci} = mdl;
    grps = [ones(size(trCFR));(2*ones(size(trAFR)))];
    mdlT = table([trCFR;trAFR],[trC;trA],[ones(size(trCFR));(2*ones(size(trAFR)))]);
    mdlT.Properties.VariableNames = {'FR','TR','GR'};
    mdlT.GR = categorical(mdlT.GR);
    fit = fitlm(mdlT,'TR~FR*GR');
    ant = anova(fit);
    pvals(ci) = ant{3,5};
    minFR = min(mdlT{:,1}); maxFR = max(mdlT{:,1});
    xFR = [minFR:0.001:maxFR]'; gFR = ones(size(xFR)); 
    xFR = [xFR;xFR]; gFR = [gFR;(2*gFR)];
    xdT = table(xFR,gFR); xdT.Properties.VariableNames = mdlT.Properties.VariableNames([1 3]); xdT.GR = categorical(xdT.GR);
    trCp = predict(fit,xdT);
    plot(xFR(gFR == 1),trCp(gFR == 1),'k');
    plot(xFR(gFR == 2),trCp(gFR == 2),'r');
%     plot(trAFR,trCp(gFR == 2),'r');
    ylim([0 max(mdlT{:,2})]); xlim([0 max(mdlT{:,1})]);
    set(gca,'FontSize',6,'FontWeight','Bold','TickDir','out','xcolor','k','ycolor','k');
    if ci == 1
        put_axes_labels(ha,{{'Avg. FR(A.U.)'},[0 0 0]},{{'Avg. # of Trans./min'},[0 0 0]});
    else
        put_axes_labels(ha,{{'Avg. FR(A.U.)'},[0 0 0]},{{''},[0 0 0]});
    end
    titletxt = sprintf('C%d',ci); ht = set_axes_top_text_no_line(gcf,ha,titletxt,[0.031 0.1 0 0]);set(ht,'FontSize',8);
    titletxt = sprintf('%s',getNumberOfAsterisks(pvals(ci))); ht = set_axes_top_text_no_line(gcf,ha,titletxt,[0.031 -0.031 0 0]);set(ht,'FontSize',8);
    format_axes_b(gca);
end
save_pdf(ff.hf,mData.pdf_folder,sprintf('Distribution_firing_rate'),600);
%% make scatter and linear regression plots between FR and RF,zMI,Rs
% Select Air ON or Air OFF phases
for piii = 1:2%length(varT)
    magfac = mData.magfac;
    ff = makeFigureRowsCols(108,[4 5 3 1],'RowsCols',[1 8],'spaceRowsCols',[0.03 -0.01],'rightUpShifts',[0.1 0.3],'widthHeightAdjustment',[10 -585]);
    MY = 10; ysp = 1; mY = 0; 
    stp = 0.4*magfac; widths = ([0.5 0.5 0.5 0.5 0.25 0.25 0.25 0.25]+0.0015)*magfac; gap = 0.175*magfac;
    adjust_axes(ff,[mY MY],stp,widths,gap,{''});
    shift_axes(ff,5:8,0.35,gap);
    ylbls = {'RF','zMI','Rs'};
    confsi = 1:4; confsi = 5:8 
    var_C = all_var_C{piii}; var_A = all_var_A{piii};
    var_CvD = var_C(:,confsi); var_AvD = var_A(:,confsi); var_CvDFR = out_C(:,confsi); var_AvDFR = out_A(:,confsi);
    for ci = 1:4
        distD = [var_CvD(:,ci) var_AvD(:,ci)]; distDFR = [var_CvDFR(:,ci) var_AvDFR(:,ci)];
        [distDo,allVals,allValsG] = plotDistributions(distD);  [distDoFR,allValsFR,allValsGFR] = plotDistributions(distDFR);
        trC = allValsG{1}; trA = allValsG{2};  trCFR = allValsGFR{1}; trAFR = allValsGFR{2};
        minBin = min(allVals);
        axes(ff.h_axes(1,ci)); hold on;
        ha = gca;
        sc = scatter(trCFR,trC,1,'.','k');% sc.MarkerEdgeAlpha = 0.5; sc.MarkerFaceAlpha = 0.5
        sc = scatter(trAFR,trA,1,'.','r');% sc.MarkerEdgeAlpha = 0.5; sc.MarkerFaceAlpha = 0.5

        opts = statset('Display','iter','TolFun',1e-10);
    %     mdl = fitnlm(trCFR,trC,modelfun,beta0,'Options',opts); mdl_C{ci} = mdl;
        grps = [ones(size(trCFR));(2*ones(size(trAFR)))];
        mdlT = table([trCFR;trAFR],[trC;trA],[ones(size(trCFR));(2*ones(size(trAFR)))]);
        mdlT.Properties.VariableNames = {'FR','TR','GR'};
        mdlT.GR = categorical(mdlT.GR);
        fit = fitlm(mdlT,'TR~FR*GR');
        ant = anova(fit);
        pvals(ci) = ant{3,5};
        minFR = min(mdlT{:,1}); maxFR = max(mdlT{:,1});
        xFR = [minFR:0.001:maxFR]'; gFR = ones(size(xFR)); 
        xFR = [xFR;xFR]; gFR = [gFR;(2*gFR)];
        xdT = table(xFR,gFR); xdT.Properties.VariableNames = mdlT.Properties.VariableNames([1 3]); xdT.GR = categorical(xdT.GR);
        trCp = predict(fit,xdT);
        plot(xFR(gFR == 1),trCp(gFR == 1),'k');
        plot(xFR(gFR == 2),trCp(gFR == 2),'r');
    %     plot(trAFR,trCp(gFR == 2),'r');
        ylim([0 max(mdlT{:,2})]); xlim([0 max(mdlT{:,1})]);
        set(gca,'FontSize',6,'FontWeight','Bold','TickDir','out','xcolor','k','ycolor','k');
        if ci == 1
            put_axes_labels(ha,{{'Avg. FR(A.U.)'},[0 0 0]},{{ylbls{piii}},[0 0 0]});
        else
            put_axes_labels(ha,{{'Avg. FR(A.U.)'},[0 0 0]},{{''},[0 0 0]});
        end
        titletxt = sprintf('C%d',ci); ht = set_axes_top_text_no_line(gcf,ha,titletxt,[0.031 0.1 0 0]);set(ht,'FontSize',8);
        titletxt = sprintf('%s',getNumberOfAsterisks(pvals(ci))); ht = set_axes_top_text_no_line(gcf,ha,titletxt,[0.031 -0.031 0 0]);set(ht,'FontSize',8);
        format_axes_b(gca);
    end
end
save_pdf(ff.hf,mData.pdf_folder,sprintf('Distribution_firing_rate'),600);

