function raster_properties

% general for all properties including responsivity, response fidelity, zMI, Rs
ntrials = 50;
si = [C1_t_D C2_t_D C3_t_D C4_t_D];
%     si = [C1_i_T C2_i_T C3_i_T C4_i_T];
Rs_C = oC.Rs(:,si); Rs_A = oA.Rs(:,si); mRs_C = oC.mR(:,si); mRs_A = oA.mR(:,si);
props_C = get_props_Rs(oC.Rs(:,si),ntrials); props_A = get_props_Rs(oA.Rs(:,si),ntrials);
% pop_var_name = {'all','vals','valsT','Nvals','good_zMI','Ngood_zMI'};
pop_var_name = {'good_zMI','good_Gauss','good_MFR'};
pop_var_name = {'all'};
% pop_var_name = {'vals'};
% pop_var_name = {'vals','good_zMI'};
sel_pop_C = cell_list_op(props_C,pop_var_name); sel_pop_A = cell_list_op(props_A,pop_var_name);
pop_var_name = {'Nvals'};
sel_pop_CNR = cell_list_op(props_C,pop_var_name); sel_pop_ANR = cell_list_op(props_A,pop_var_name);


pop_var_name = {'vals'};
rsel_pop_C = cell_list_op(props_C,pop_var_name); rsel_pop_A = cell_list_op(props_A,pop_var_name);

cell_types = {'C1','C2','C3','C4'};
bins = [0 50 100 150]; 
bins = [0 150];
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
mean_var_C = [];mean_var_A = [];

params = {'perc','N_Resp_Trials','zMI','rs','nan_zMI','nan_rs','HaFD','HiFD','PWs','centers','peak_locations','mean_FR','MFR'};
varT = 3;%:length(params)
for pii = varT
    if pii == 1
        mean_var_C = exec_fun_on_cell_mat(sel_pop_C,'percent'); mean_var_A = exec_fun_on_cell_mat(sel_pop_A,'percent');
    else
        eval(sprintf('var_C = get_vals(props_C.%s,sel_pop_C);',params{pii})); eval(sprintf('var_A = get_vals(props_A.%s,sel_pop_A);',params{pii}));
        if pii == 5 || pii == 6
            mean_var_C = exec_fun_on_cell_mat(sel_pop_C,'percent'); mean_var_A = exec_fun_on_cell_mat(sel_pop_A,'percent'); 
        else
            mean_var_C = exec_fun_on_cell_mat(var_C,'nanmean'); mean_var_A = exec_fun_on_cell_mat(var_A,'nanmean'); 
        end
    end
end

all_vars_C = []; all_vars_A = [];
for ii = 1:4
    t_all_varsC = []; t_all_varsA = [];
    for ani = 1:5
        t_all_varsC = [t_all_varsC;var_C{ani,ii}];
        t_all_varsA = [t_all_varsA;var_A{ani,ii}];
    end
    all_vars_C(:,ii) = t_all_varsC; all_vars_A(:,ii) = t_all_varsA;
end

% all_vars_C = array2table(all_vars_C); all_vars_A = array2table(all_vars_A);
% all_vars_C.Properties.VariableNames = {'C1','C2','C3','C4'}; all_vars_A.Properties.VariableNames = {'C1','C2','C3','C4'};
% all_vars_C(sum(isnan(all_vars_C), 2) == 1, :) = [];

[within,dvn,xlabels] = make_within_table({'Cond'},[4]);
[dataT,dataTNC] = make_between_table({all_vars_C;all_vars_A},dvn);
ra = RMA(dataT,within,{0.05,{'bonferroni'}});
ra.ranova
print_for_manuscript(ra)
%%
for c1i = 1:4
    for c2i = 1:4
        scni = [c1i c2i];
        [within,dvn,xlabels] = make_within_table({'Cond'},[2]);
        dataT = make_between_table({all_vars_C(:,scni);all_vars_A(:,scni)},dvn);
        ra = RMA(dataT,within,{0.05,{'bonferroni'}});
%         ra.ranova
%         print_for_manuscript(ra)
        groupps(c1i,c2i) = ra.ranova{2,ra.selected_pval_col};
        condps(c1i,c2i) = ra.ranova{4,ra.selected_pval_col};
        groupcondps(c1i,c2i) = ra.ranova{5,ra.selected_pval_col};
    end
end

%% NLME mixed effects modeling I am trying to fit sinusoidal model
num = table2array(dataTNC); naninds = sum(isnan(num),2)>0; num(naninds,:) = [];
adata = make_text_file_for_nlme_R(num,'./AD_Paper_Revamp/nlme_zMI.txt');
group = categorical(adata(:,4));
dv = dummyvar(group);
X = [adata(:,2) dv]; y = adata(:,3); X(isnan(y),:) = []; y(isnan(y)) = [];
modelfun = @(b,x)X(:,2).*b(1).*sin(2.*pi.*b(2).*X(:,1)+b(3)) + X(:,3).*b(4).*sin(2.*pi.*b(5).*X(:,1)+b(6));
beta0 = [1 1 1 1 1 1]; % for each parameter, define values for both the groups 
mdl = fitnlm(X,y,modelfun,beta0);

%%
X = adata(:,2); y = adata(:,3); group = adata(:,4);
% yfit = modelfun([60 0.1],,VFUN)
model = @(b,x)b(1).*sin(2.*pi.*b(2).*x+b(3));
beta0 = [1 1 1]; fegd(:,:,1) = beta0; fegd(:,:,2) = beta0;
options = statset('nlmefit');
clc
group = categorical(group);
options = statset(options,'MaxIter',350,'TolX',1e-8,'Display','final');
[beta,PSI,stats,B] = nlmefit(X,y,group,[],model,beta0,'Options',options,'RefineBeta0','On');


%%
combs = {[1 2],[2 3],[3 4]}; dp = [];
for ani = 1:5
    for ii = 1:length(combs)
        tcomb = combs{ii};
        tC = cell_list_op(rsel_pop_C(ani,combs{ii}),[],'and',1); tA = cell_list_op(rsel_pop_A(ani,combs{ii}),[],'and',1);
        tvC1 = var_C{ani,tcomb(1)}; tvC2 = var_C{ani,tcomb(2)}; tvA1 = var_A{ani,tcomb(1)}; tvA2 = var_A{ani,tcomb(2)};
        p1C = tvC1(tC{1}); p2C = tvC2(tC{1}); p1A = tvA1(tA{1}); p2A = tvA2(tA{1});
        dpC{ani,ii} = p2C-p1C; dpA{ani,ii} = p2A-p1A;
        numC(ani,ii) = length(p1C); numA(ani,ii) = length(p1A);
        totC(ani,ii) = length(tC{1}); totA(ani,ii) = length(tA{1});
    end
end

varC = 100 * numC./totC;
varA = 100 * numA./totA;
varC = exec_fun_on_cell_mat(dpC,'mean');
varA = exec_fun_on_cell_mat(dpA,'mean');
[within,dvn,xlabels] = make_within_table({'Cond'},[3]);
dataT = make_between_table({varC;varA},dvn);
ra = RMA(dataT,within,{0.05,{'hsd'}});
ra.ranova
print_for_manuscript(ra)

%%

typeCorr = {'Spatial Correlation',{'Population Vector','Correlation'},'\Delta FR Score'};
FF = {'SP','PV','RR'};
typeCorr = {'Temporal Correlation',{'Population Vector','Correlation'},'\Delta FR Score'};
FF = {'TP','PV','RR'};
varT = 1;
xlabels = {titletxt,titletxt,titletxt};
ylabels = xlabels;
incrs = [0.01,0.01,0.01];
mYs = [-0.1,0,0];
MYs = [0.25,0.2,1.2];
ysps = [0.04,0.04,0.1];

incrs = [0.01,0.01,0.01];
mYs = [-1.5,0,0];
MYs = [2.5,0.2,1.2];
ysps = [0.03,0.04,0.1];

magfac = mData.magfac;
ff = makeFigureRowsCols(108,[4 5 6.9 1.45],'RowsCols',[1 8],'spaceRowsCols',[0.03 -0.01],'rightUpShifts',[0.1 0.3],'widthHeightAdjustment',[10 -585]);
% MY = 8; ysp = 1; mY = 0; 
stp = 0.4*magfac; widths = ([0.5 0.5 0.5 0.5 0.25 0.25 0.25 0.25]+0.23)*magfac; gap = 0.175*magfac;
adjust_axes(ff,[mY MY],stp,widths,gap,{'Cell #'});

shift_axes(ff,5:8,0.35,gap);
var_CvD = dpC; var_AvD = dpA;
for ci = 1:3
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
    [ha,hb,hca] = plotDistributions(allValsG,'colors',tcolors,'maxY',100,'min',minBin,'incr',incr,'max',maxBin,'do_mean','No');
%     [ha,hb,hca] = plotDistributions(distDo,'colors',tcolors,'maxY',100,'min',minBin,'incr',incr,'max',maxBin,'do_mean','Yes');
    set(gca,'FontSize',6,'FontWeight','Bold','TickDir','out','xcolor','k','ycolor','k');
%     changePosition(gca,[0.129 0.15 -0.09 -0.13]);
    ylim([0 100]); xlim([minBin maxBin]);
    if ci == 1
        if ~isempty(strfind(xlabels{1},'PV'))
            put_axes_labels(ha,{xlabels{varT},[0 0 0]},{{'Bins (%)'},[0 0 0]});
        else
            put_axes_labels(ha,{xlabels{varT},[0 0 0]},{{'Cells (%)'},[0 0 0]});
        end
    else
%         put_axes_labels(ha,{'Firing Rate (A.U.)',[0 0 0]},{{'Neurons (%)'},[0 0 0]});
        put_axes_labels(ha,{xlabels{varT},[0 0 0]},{{''},[0 0 0]});
    end
    format_axes_b(ha);
    [ks2.h,ks2.p,ks2.ks2stat] = kstest2(allValsG{1},allValsG{2}); ks2.DF1 = length(allValsG{1}); ks2.DF2 = length(allValsG{2});
    print_for_manuscript(ks2,'KS2');
    
    ht = set_axes_top_text_no_line(gcf,ha,'KS-Test',[0.0 -0.01 0 0]);set(ht,'FontSize',7);
    titletxt = sprintf('%s',getNumberOfAsterisks(ks2.p));
    ht = set_axes_top_text_no_line(gcf,ha,titletxt,[0.061 -0.01 0 0]);set(ht,'FontSize',9);
    titletxt = sprintf('C%d%d',ci,ci+1);
    ht = set_axes_top_text_no_line(gcf,ha,titletxt,[0.031 0.1 0 0]);set(ht,'FontSize',8);
    if varT == 1 || varT == 2
        titletxt = sprintf('n = %d,',length(allValsG{1}));
        ht = set_axes_top_text_no_line(gcf,ha,titletxt,[0.04 -0.3 0 0]);set(ht,'FontSize',7,'Color','k');
        titletxt = sprintf('%d',length(allValsG{2}));
        ht = set_axes_top_text_no_line(gcf,ha,titletxt,[0.06 -0.4 0 0]);set(ht,'FontSize',7,'Color','r');
    end
    if varT == 3
        titletxt = sprintf('n = %d,',length(allValsG{1}));
        ht = set_axes_top_text_no_line(gcf,ha,titletxt,[0.01 -0.3 0 0]);set(ht,'FontSize',7,'Color','k');
        titletxt = sprintf('%d',length(allValsG{2}));
        ht = set_axes_top_text_no_line(gcf,ha,titletxt,[0.03 -0.4 0 0]);set(ht,'FontSize',7,'Color','r');
    end
end
for ci = 1:3
    distD = [var_CvD(:,ci) var_AvD(:,ci)];
    tcolors = {'k','r'};
    [distDo,allVals,allValsG] = plotDistributions(distD);
    axes(ff.h_axes(1,ci+4)); hold on;
    mVar = []; semVar = [];
    mVar(1,1) = nanmean(allValsG{1}); mVar(2,1) = nanmean(allValsG{2});
    semVar(1,1) = nanstd(allValsG{1})/sqrt(length(allValsG{1})); semVar(2,1) = nanstd(allValsG{2})/sqrt(length(allValsG{2}));
    mY = mYs(varT); MY = MYs(varT); ysp = ysps(varT);%max([mVar+semVar]);
    tcolors = {'k','r'};%{colors{1};colors{2};colors{3};colors{4};colors{1};colors{2};colors{3};colors{4}};
    xdata = make_xdata([2],[1 2]);   combs = [1 2];
    
    [t2.h,t2.p,t2.ci,t2.tstat] = ttest2(allValsG{1},allValsG{2}); t2.cd = computeCohen_d(allValsG{1},allValsG{2});
    print_for_manuscript(t2,'t2');
    
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[t2.h t2.p],'colors',tcolors,'sigColor','k',...
        'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',9,'barWidth',0.5,'sigLinesStartYFactor',0.05);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[mY MY]); xticks = xdata; 
    if ci == 1
        hy = ylabel(ylabels{varT}); changePosition(hy,[-0.5 0 0]);
        set(gca,'ytick',[mY MY],'yticklabels',{sprintf('%.1f',mY),sprintf('%.1f',MY)});
    end
    xticklabels = {'C-TG','A-TG'};set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45);
    make_bars_hollow(hbs(1:end));
%     put_axes_labels(gca,{'',[]},{ylabeltxt,[]});
    format_axes_b(gca);
    titletxt = sprintf('C%d%d',ci,ci+1);
    ht = set_axes_top_text_no_line(gcf,gca,titletxt,[0.021 0.1 0 0]);set(ht,'FontSize',8);
    ht = set_axes_top_text_no_line(gcf,gca,'t-Test',[0.0 -0.01 0 0]);set(ht,'FontSize',7);
end
delete(ff.h_axes(1,[4 8]))
save_pdf(ff.hf,mData.pdf_folder,sprintf('Distribution_firing_rate'),600);
