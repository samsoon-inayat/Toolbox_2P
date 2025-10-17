%% first find the variables time_cells, distance_cells, and speed_cells in the final_time_dist_bin...m file
variable_combs = {'FR_time','FR_dist','FR_speed'};
metric = 'MI';
metric = 'PC';
clear time_cells distance_cells speed_cells
clear t_cells_T d_cells_T s_cells_T td_cells_T ds_cells_T ts_cells_T tds_cells_T ntds_cells_T
clear t_cells_I d_cells_I s_cells_I td_cells_I ds_cells_I ts_cells_I tds_cells_I ntds_cells_I
all_cells = {};
si = [Ar_t_T Ar_i_T Ar_t_D Ar_i_D ArL_t_T ArL_i_T ArL_t_D ArL_i_D Ars_t_T Ars_i_T Ars_t_D Ars_i_D]; propsPL = get_props_Rs(o.Rs(:,si),30);
si = [Ar_t_T Ar_i_T ArL_t_T ArL_i_T Ars_t_T Ars_i_T]; si_cn_ap = [[1 1 2 2 3 3];[1 2 1 2 1 2]];
sib = [Ab_t_T Ab_i_T Abs_t_T Abs_i_T]; propsB = get_props_Rs(o.Rs(:,sib),30);
% si = [Ar_t_D Ar_i_D ArL_t_D ArL_i_D Ars_t_D Ars_i_D]; 
props = get_props_Rs(o.Rs(:,si),30);
propsT = get_props_new(outT,met_valsT,props,si_cn_ap);
propsD = get_props_new(outD,met_valsD,props,si_cn_ap);
% propsTDM = get_props_new_matched({outT;outD},{met_valsT;met_valsD},props,si_cn_ap);


siB = [Ab_t_T Ab_i_T Abs_t_T Abs_i_T]; si_cn_apB = [[1 1 2 2];[1 2 1 2]];
propsTB = get_props_newB(outTB,met_valsTB,propsB,si_cn_apB);
% [SinT,MixT,AllT] = get_the_pops(propsT,propsD);
%%
cell_popsPC = [];%[propsT.newPC.cells_time propsD.newPC.cells_dist propsT.newPC.cells_speed];
cell_popsMI = [];%[propsT.newMI.cells_time propsD.newMI.cells_dist propsT.newMI.cells_speed];
for an = 1:5
    cni = 1;
    for cc = 1:6
        cell_popsPC{an,cni} = propsT.newPC.cells_time{an,cc};
        cell_popsMI{an,cni} = propsT.newMI.cells_time{an,cc};
        cni = cni + 1;
        cell_popsPC{an,cni} = propsD.newPC.cells_dist{an,cc};
        cell_popsMI{an,cni} = propsD.newMI.cells_dist{an,cc};
        cni = cni + 1;
        cell_popsPC{an,cni} = propsT.newPC.cells_speed{an,cc};
        cell_popsMI{an,cni} = propsT.newMI.cells_speed{an,cc};
        cni = cni + 1;
    end
end

cell_popsPCB = [];%[propsT.newPC.cells_time propsD.newPC.cells_dist propsT.newPC.cells_speed];
cell_popsMIB = [];%[propsT.newMI.cells_time propsD.newMI.cells_dist propsT.newMI.cells_speed];
for an = 1:5
    cni = 1;
    for cc = 1:4
        cell_popsPCB{an,cni} = propsTB.newPC.cells_time{an,cc};
        cell_popsMIB{an,cni} = propsTB.newMI.cells_time{an,cc};
        cni = cni + 1;
        cell_popsPCB{an,cni} = propsTB.newPC.cells_speed{an,cc};
        cell_popsMIB{an,cni} = propsTB.newMI.cells_speed{an,cc};
        cni = cni + 1;
    end
end

%%
clc
eq = [];
for an = 1:5
    for cc = 1:18
        eq(an,cc) = isequal(pop.X{an,cc},cell_popsPC{an,cc});
    end
end
eq



%%
MVT = met_valsT; MVD = met_valsD;
for an = 1:5
    for cni = 1:size(si_cn_ap,2)
        cn = si_cn_ap(1,cni); ap = si_cn_ap(2,cni);
        idx = 3; pvalsPC = [MVT{1}{an,cn,ap}.PC(:,idx) MVD{2}{an,cn,ap}.PC(:,idx) MVT{3}{an,cn,ap}.PC(:,idx)];
        idx = 3; pvalsMI = [MVT{1}{an,cn,ap}.MI(:,idx) MVD{2}{an,cn,ap}.MI(:,idx) MVT{3}{an,cn,ap}.MI(:,idx)];
        idvPC{an,cni} = pvalsPC < 0.05;
        idvMI{an,cni} = pvalsMI < 0.05;
        pvalsPC_a{an,cni} = pvalsPC;
        pvalsMI_a{an,cni} = pvalsMI;
        idx = 2; zvalsPC = [MVT{1}{an,cn,ap}.PC(:,idx) MVD{2}{an,cn,ap}.PC(:,idx) MVT{3}{an,cn,ap}.PC(:,idx)];
        idx = 2; zvalsMI = [MVT{1}{an,cn,ap}.MI(:,idx) MVD{2}{an,cn,ap}.MI(:,idx) MVT{3}{an,cn,ap}.MI(:,idx)];
        zvalsPC_a{an,cni} = zvalsPC;
        zvalsMI_a{an,cni} = zvalsMI;
    end
end

pvals = pvalsPC_a; zvals = zvalsPC_a;
disp('Done')


%%
clc
% Use either idv = idvPC or idvMI
idv = idvMI;  A = size(idv,1);  C = 6;
delta = nan(A,1);

for a = 1:A
    [~,~,~,~,P] = jaccard_animal(idv,a);   % agreement matrix, diag=NaN
    same = false(C); same(1:2:end,1:2:end)=true; same(2:2:end,2:2:end)=true; % A1/A1 & A2/A2
    offdiag = true(C); offdiag(1:C+1:end)=false;

    m_same = mean(P(same & offdiag),'omitnan');
    m_diff = mean(P(~same & offdiag),'omitnan');
    delta(a) = m_same - m_diff;
end

fprintf('Δ per animal: '); fprintf('%.3f ', delta); fprintf('\n');
fprintf('Group: mean Δ = %.3f ± %.3f SEM\n', mean(delta), std(delta)/sqrt(A));

[pW,~,~] = signrank(delta,0,'tail','right');   % expect Δ>0
fprintf('Wilcoxon signed-rank (Δ>0): p = %.4g\n', pW);

%% Set input (you already have these in workspace)
% idv = idvPC;              % or: idv = idvMI;

txl_all6 = []; cnti = 1;
for cn = 1:3
    for ap = 1:2
        txl_all6{cnti} = sprintf('C%d-A%d',cn+2,ap);
        cnti = cnti + 1;
    end
end
ct = {'T','D','S'};
txl_all18 = []; cnti = 1;
for cn = 1:3
    for ap = 1:2
        for cti = 1:3
            txl_all18{cnti} = sprintf('C%d-A%d-%s',cn+2,ap,ct{cti});
            cnti = cnti + 1;
        end
    end
end

ct = {'T','S'};
txl_all8B = []; cnti = 1;
cnums = [2 7];
for cn = 1:2
    for ap = 1:2
        for cti = 1:2
            txl_all8B{cnti} = sprintf('C%d-A%d-%s',cnums(cn),ap,ct{cti});
            cnti = cnti + 1;
        end
    end
end


%% --- Pick one animal for detailed view
idv = idvPC;
txl = txl_all6;
G = jaccard_group(idv);
a = 3;
[D_cell, D_mean, D_sem, N_pair, P_same] = jaccard_animal(idv, a);
adj_idx = sub2ind([6 6], 1:5, 2:6);
mean_adj = mean(D_mean(adj_idx), 'omitnan');
fprintf('Animal %d: mean adjacent identity distance = %.3f\n', a, mean_adj);

ff = makeFigureRowsCols(107,[5 3 3.65 2.25],'RowsCols',[2 4],'spaceRowsCols',[0.2 0.0451],'rightUpShifts',[0.01 0.119],'widthHeightAdjustment',[-10 -200]);
    set(gcf,'color','w');     ylims = [0 1];
    stp = 0.2575; widths = 0.35*([2 2 2 2 0.4 0.4]+1.05); gap = -0.30251572191; adjust_axes(ff,ylims,stp,widths,gap,{'Euclidean Distance'});
    gap1 = 0.27; widths1 = 1*0.75; height1 = widths1 - 0.05;

cmaplims = [0 0.35];
axes(ff.h_axes(1,1)); imt = imagesc(P_same,cmaplims); axis square; colormap(parula); %colorbar
format_axes(gca); set(gca,'xtick',1:length(txl),'ytick',1:length(txl_all6),'xticklabels',txl_all6,'yticklabels',txl_all6,'Ydir','normal'); xtickangle(45);
ytickangle(45);hyl = ylabel({'Animal 3','N=315 Cells'}); changePosition(hyl,[0.5 0 0])
imAlpha=ones(size(P_same));imAlpha(isnan(P_same)) = 0;imt.AlphaData = imAlpha;set(gca,'Color',[0.75 0.75 0.75])
axes(ff.h_axes(1,2)); imt = imagesc(D_sem,cmaplims); axis square; colormap(parula); %colorbar
format_axes(gca); set(gca,'xtick',1:length(txl),'ytick',1:length(txl_all6),'xticklabels',txl_all6,'yticklabels','','Ydir','normal'); xtickangle(45);
ytickangle(45); 
imAlpha=ones(size(P_same));imAlpha(isnan(P_same)) = 0;imt.AlphaData = imAlpha;set(gca,'Color',[0.75 0.75 0.75])

% title('SEM over cells');

axes(ff.h_axes(2,1)); imt = imagesc(G.P_group_mean,cmaplims); axis square; colormap(parula); %colorbar
format_axes(gca); set(gca,'xtick',1:length(txl),'ytick',1:length(txl_all6),'xticklabels',txl_all6,'yticklabels',txl_all6,'Ydir','normal'); xtickangle(45);
ytickangle(45);hyl = ylabel({'Group','N = 5 Animals'}); changePosition(hyl,[0.5 0 0])
% title('Agreement = 1 - distance');
imAlpha=ones(size(P_same));imAlpha(isnan(P_same)) = 0;imt.AlphaData = imAlpha;set(gca,'Color',[0.75 0.75 0.75])

axes(ff.h_axes(2,2)); imt = imagesc(G.D_group_sem,cmaplims); axis square; colormap(parula); %colorbar
format_axes(gca); set(gca,'xtick',1:length(txl),'ytick',1:length(txl_all6),'xticklabels',txl_all6,'yticklabels','','Ydir','normal'); xtickangle(45);
ytickangle(45);
imAlpha=ones(size(P_same));imAlpha(isnan(P_same)) = 0;imt.AlphaData = imAlpha;set(gca,'Color',[0.75 0.75 0.75])

% title('SEM over cells');


idv = idvMI;
G = jaccard_group(idv);
a = 3;
[D_cell, D_mean, D_sem, N_pair, P_same] = jaccard_animal(idv, a);
adj_idx = sub2ind([6 6], 1:5, 2:6);
mean_adj = mean(D_mean(adj_idx), 'omitnan');
fprintf('Animal %d: mean adjacent identity distance = %.3f\n', a, mean_adj);

cmaplims = [0 0.35];
axes(ff.h_axes(1,3)); imt = imagesc(P_same,cmaplims); axis square; colormap(parula); %colorbar
format_axes(gca); set(gca,'xtick',1:length(txl),'ytick',1:length(txl_all6),'xticklabels',txl_all6,'yticklabels','','Ydir','normal'); xtickangle(45);
ytickangle(45);
imAlpha=ones(size(P_same));imAlpha(isnan(P_same)) = 0;imt.AlphaData = imAlpha;set(gca,'Color',[0.75 0.75 0.75])

% title('Agreement = 1 - distance');
axes(ff.h_axes(1,4)); imt = imagesc(D_sem,cmaplims); axis square; colormap(parula); %colorbar
format_axes(gca); set(gca,'xtick',1:length(txl),'ytick',1:length(txl_all6),'xticklabels',txl_all6,'yticklabels','','Ydir','normal'); xtickangle(45);
ytickangle(45);
imAlpha=ones(size(P_same));imAlpha(isnan(P_same)) = 0;imt.AlphaData = imAlpha;set(gca,'Color',[0.75 0.75 0.75])

[hc,ha] = putColorBar(ff.h_axes(1,4),[-0.025 0.061 -0.21 -0.15],cmaplims,6,'eastoutside',[0.1 0.075 0.1 0.075]);
% title('SEM over cells')5

axes(ff.h_axes(2,3)); imt = imagesc(G.P_group_mean,cmaplims); axis square; colormap(parula); %colorbar
format_axes(gca); set(gca,'xtick',1:length(txl),'ytick',1:length(txl_all6),'xticklabels',txl_all6,'yticklabels','','Ydir','normal'); xtickangle(45);
ytickangle(45);
imAlpha=ones(size(P_same));imAlpha(isnan(P_same)) = 0;imt.AlphaData = imAlpha;set(gca,'Color',[0.75 0.75 0.75])

% title('Agreement = 1 - distance');
axes(ff.h_axes(2,4)); imt = imagesc(G.D_group_sem,cmaplims); axis square; colormap(parula); %colorbar
format_axes(gca); set(gca,'xtick',1:length(txl),'ytick',1:length(txl_all6),'xticklabels',txl_all6,'yticklabels','','Ydir','normal'); xtickangle(45);
ytickangle(45);
imAlpha=ones(size(P_same));imAlpha(isnan(P_same)) = 0;imt.AlphaData = imAlpha;set(gca,'Color',[0.75 0.75 0.75])
[hc,ha] = putColorBar(ff.h_axes(2,4),[-0.025 0.061 -0.21 -0.15],cmaplims,6,'eastoutside',[0.1 0.075 0.1 0.075]);
% title('SEM over cells');

for rr = 1:2
    for cc = 1:4
        axes(ff.h_axes(rr,cc));
        if ismember(cc,[1 3])
            set_axes_top_text_no_line(ff.hf,ff.h_axes(rr,cc),'Mean',[0.1 0.04 0 0]);
        else
            set_axes_top_text_no_line(ff.hf,ff.h_axes(rr,cc),'SEM',[0.1 0.04 0 0]);
        end
    end
end

save_pdf(ff.hf,mData.pdf_folder,sprintf('AG_Map.pdf'),600);

%% --- Group Jaccard stats
figure('Color','w');
subplot(1,2,1); imagesc(G.D_group_mean,[0 1]); axis square; colormap(parula); colorbar
title('Group Jaccard distance (mean over animals)');
subplot(1,2,2); imagesc(G.P_group_mean,[0 1]); axis square; colormap(parula); colorbar
title('Group agreement (=1 - distance)');

%% --- Stability percentages (Stable / Dynamic / Middle per label)
P = stability_percentages(idv);
labels = {'Time','Distance','Speed'};
fprintf('\nPer-animal %%Stable / %%Dynamic / %%Middle (denominators = cells that ever had the label)\n');
for a = 1:size(idv,1)
    fprintf('Animal %d:\n', a);
    for b = 1:3
        fprintf('  %-8s  S=%.1f%%  D=%.1f%%  M=%.1f%%   (n=%d)\n', ...
            labels{b}, P.stable(a,b), P.dynamic(a,b), P.middle(a,b), P.denom(a,b));
    end
end

fprintf('\nGroup means ± SEM (%%):\n');
for b = 1:3
    fprintf('  %-8s  Stable  %5.1f ± %4.1f   Dynamic %5.1f ± %4.1f   Middle %5.1f ± %4.1f\n', ...
        labels{b}, ...
        P.mean_sem.stable_mean(b),  P.mean_sem.stable_sem(b), ...
        P.mean_sem.dynamic_mean(b), P.mean_sem.dynamic_sem(b), ...
        P.mean_sem.middle_mean(b),  P.mean_sem.middle_sem(b));
end
%%
    [OIo,mOIo,semOIo,OI_mat,p_vals,h_vals,all_CI,mCI,semCI,all_CI_mat,uni] = get_overlap_index(cell_popsMI,0.5,0.05);
    % mOI = mCI; semOI = semCI;
    % mCI = mOIo; semCI = semOIo;

    mUni = nanmean(uni,3); mUni1 = tril(mUni,-1) + tril(mUni,-1)'; mUni2 = triu(mUni,1) + triu(mUni,1)'; mmUni = min(mUni(:)); MmUni = max(mUni(:));
    semUni = nanstd(uni,[],3)/sqrt(5); semUni1 = tril(semUni,-1) + tril(semUni,-1)'; semUni2 = triu(semUni,1) + triu(semUni,1)'; msemUni = min(semUni(:)); MsemUni = max(semUni(:));
%     figure(1000);clf;subplot 131;imagesc(mUni,[mmUni MmUni]);set(gca,'YDir','normal');subplot 132;imagesc(mUni1,[mmUni MmUni]);set(gca,'YDir','normal');subplot 133;imagesc(mUni2,[mmUni MmUni]);set(gca,'YDir','normal');
    
   ff = makeFigureRowsCols(107,[5 3 6.9 2.5],'RowsCols',[1 2],'spaceRowsCols',[0.2 0.2],'rightUpShifts',[0.12 0.19],'widthHeightAdjustment',[-10 -300]);
    set(gcf,'color','w');     ylims = [0 1];
    stp = 0.57; widths = 0.5*([2 2 2 2 0.4 0.4]+1.5); gap = 0.99*1.75; adjust_axes(ff,ylims,stp,widths,gap,{'Euclidean Distance'});
    gap1 = 0.32; widths1 = 1*0.75; height1 = widths1 - 0.05;
    mOIo = mOIo(ord,ord); txl = txl_all18(ord);
    mats = {mOIo,mUni}; semmats = {semOIo,semUni};
    titles = {'Mean','Complementation'};
    for ii = 1%:length(mats)
        mOI = mats{ii}; semOI = semmats{ii};
        sz = size(mOI,1);        oM = ones(size(mOI));
        mask1 = (triu(oM,0) & tril(oM,0)); mOI(mask1==1) = NaN; semOI(mask1 == 1) = NaN;
        maxI = max([mOI(:);semOI(:)]);            minI = min([mOI(:);semOI(:)]);

        mask = tril(NaN(size(mOI)),0); mask(mask==0) = 1; 

        imAlpha=ones(size(mOI));    %imAlpha(isnan(mask))=0.25; 
        imAlpha(mask1 == 1) = 0;
        axes(ff.h_axes(1,ii));
        im1 = imagesc(mOI,[minI,maxI]);    im1.AlphaData = imAlpha;
        set(gca,'color',0.5*[1 1 1]);    colormap parula;    %axis equal
        format_axes(gca);
        set_axes_limits(gca,[0.5 sz+0.5],[0.5 sz+0.5]);
        set(gca,'xtick',1:length(txl),'ytick',1:length(txl),'xticklabels',txl,'yticklabels',txl,'Ydir','reverse'); xtickangle(45);
        for rr = 1:size(mCI,1)
            for cc = 1:size(mCI,1)
                if rr == cc
                    % text(rr-0.25,cc,sprintf('%.0f',mCI(rr,cc)),'FontSize',5);
                end
            end
        end
    %     plot([(10.5+((ii-1)*10)) (10.5+((ii-1)*10))],[0 150.5],'w','linewidth',0.1); 
    %     plot([0 150.5],[(10.5+((ii-1)*10)) (10.5+((ii-1)*10))],'w','linewidth',0.1); 
        set(gca,'Ydir','normal');ytickangle(15);        box on
        format_axes(gca);
        [hc,hca] = putColorBar(gca,[0.09 0.07 -0.11 -0.15],{sprintf('%.1f',minI),sprintf('%.1f',maxI)},6,'eastoutside',[0.07 0.05 0.07 0.05]);
        colormap parula
        htit = title('Mean'); changePosition(htit,[0 0.057 0]); set(htit,'FontSize',6,'FontWeight','normal');
        % axes(ff.h_axes(1,ii));
        % ht = set_axes_top_text_no_line(gcf,gca,titles{ii},[0 0 0 0]);set(ht,'FontSize',6,'FontWeight','normal')

%         save_pdf(hf,mData.pdf_folder,sprintf('OI_Map_%d_mean_%d.pdf',ntrials,sh),600);
        hai = ff.h_axes(1,ii);
        pos = get(hai,'Position'); units = get(hai,'Units');
        ha = axes('Position',pos,'Visible','on','Units',units); poshca =  get(hca,'Position');
        set(ha,'Position',[hai.Position(1)+hai.Position(3)+gap1 (hai.Position(2)+hai.Position(4)-height1) widths1 height1]);
%         axes(ff.h_axes(1,jj(ii)+ii+1));
        im1 = imagesc(semOI,[minI,maxI]);    im1.AlphaData = imAlpha;
        set(gca,'color',0.5*[1 1 1]);    colormap parula;    %axis equal
        format_axes(gca);
        set_axes_limits(gca,[0.5 sz+0.5],[0.5 sz+0.5]);
        set(gca,'xtick',1:length(txl),'ytick',1:length(txl),'xticklabels',[],'yticklabels',[],'Ydir','reverse'); xtickangle(75);
        changePosition(gca,[0.0 0 -0.05 0]);
        set(gca,'Ydir','normal');
        htit = title('SEM'); changePosition(htit,[0 0.07 0]); set(htit,'FontSize',6,'FontWeight','normal');
        box on
        format_axes(gca);
        colormap parula
    end
    save_pdf(ff.hf,mData.pdf_folder,sprintf('OI_Map.pdf'),600);
   %% for clustering
   ff = makeFigureRowsCols(107,[5 3 6.9 2.5],'RowsCols',[2 3],'spaceRowsCols',[0.3 -0.02],'rightUpShifts',[1 0.18],'widthHeightAdjustment',[-10 -300]);
    set(gcf,'color','w');     ylims = [0 1];
    stp = 0.25; widths = [2 2 2 0.4 0.4 0.4]+0.75; gap = 0.18; adjust_axes(ff,ylims,stp,widths,gap,{'Euclidean Distance'});
    stp = 0.35; widths = 0.5*([2 2 2 0.4 0.4 0.4]+1.95); gap = 0.2; adjust_axes(ff,ylims,stp,widths,gap,{'Euclidean Distance'});

       ff = makeFigureRowsCols(107,[5 3 6.9 1.25],'RowsCols',[1 1],'spaceRowsCols',[0.3 -0.02],'rightUpShifts',[1 0.18],'widthHeightAdjustment',[-10 -300]);
    set(gcf,'color','w');     ylims = [0 1];
    stp = 0.25; widths = [2 2 2 0.4 0.4 0.4]+0.75; gap = 0.18; adjust_axes(ff,ylims,stp,widths,gap,{'Euclidean Distance'});
    stp = 0.35; widths = 0.5*([2 2 2 0.4 0.4 0.4]+1.95); gap = 0.2; adjust_axes(ff,ylims,stp,widths,gap,{'Euclidean Distance'});
    %
    hms = {1-mOIo,mUni1,mUni2};
    ahc_col_th = 0.7;
    
    hms1 = {1-mOIo,mUni1,mUni2};
    for hi = 1%:length(hms1)
        mOI1 = hms1{hi}; mOI1(mask1==1) = 0; mOI1(isnan(mOI1))=0; Di = squareform(mOI1,'tovector');%pdist(mOI1,@naneucdist); 
        % Di = pdist(mOI1,@naneucdist);
        tree = linkage(Di,'average'); [c,d] = cophenet(tree,Di); r = corr(Di',d','type','spearman');    leafOrder = optimalleaforder(tree,Di);
        hf = figure(100000000); %     leafOrder1 = leafOrder([1:3 10:12 4:9]);
    %     leafOrder1 = circshift(leafOrder,3);
        figure(hf);clf
        [H,T,TC] = dendrogram(tree,'Orientation','top','ColorThreshold',ahc_col_th*max(tree(:,3)),'Reorder',leafOrder);ylims = ylim; xlims = xlim;    set(gca,'xtick',1:length(txl),'ytick',[]);
        set(gcf,'units','inches'); set(gcf,'Position',[5 2 0.9 0.5]); close(hf);
        axes(ff.h_axes(1,hi));
        [H,T,TC] = dendrogram(tree,'Orientation','top','ColorThreshold',ahc_col_th*max(tree(:,3)),'Reorder',leafOrder);
        set(H,'linewidth',0.5); set(gca,'xtick',1:length(txl),'xticklabels',txl(leafOrder));xtickangle(45); format_axes(gca);
        if hi == 1
            hx = ylabel({'Euc. Dist.'});changePosition(hx,[0 0 0]);
        end
    %     xlim([xlims(1)+0.5 xlims(2)-0.5]);
%         changePosition(gca,[0.0 0.0 0.07 0.05]); text(0.5,ylims(2)+0,sprintf('CC = %.2f',c),'FontSize',6);
        hcct = set_axes_top_text_no_line(ff.hf,gca,sprintf('CC = %.2f',c),[0 0.05 0 0]); set(hcct,'FontWeight','Normal')
    end

    save_pdf(ff.hf,mData.pdf_folder,sprintf('OI_Map_cluster.pdf'),600);


    %%
    % === Set metric ===
txl = [txl_all8B txl_all18]
idv = idvPC;   % later: idv = idvMI;

% Optional nice labels for the 6 cases:
case_labels = {'C3-A1','C3-A2','C4-A1','C4-A2','C5-A1','C5-A2'};

% Build 18 singular-only populations and compute agreement + clustering
[pop, labels18] = build_pop18_singular(idv, case_labels);
labels18 = txl;
pop.M = 26;
pop.X = [cell_popsPCB cell_popsPC];
pop.X = [cell_popsMIB cell_popsMI];
S = pop18_agreement_cluster(pop);
ord = S.leafOrder;
% ord = 1:18;

% Plot group mean agreement (diag masked), with NaNs as light gray
A18 = S.J_mean(ord,ord);
hf = figure(1000);
set(hf,'Color','w','Position',[100 100 900 360]);
subplot(1,2,1);
h = imagesc(A18, [0 0.1]); axis square; colormap(parula); colorbar
set(h,'AlphaData',~isnan(A18)); set(gca,'Color',[0.9 0.9 0.9]);
set(gca,'XTick',1:length(ord),'XTickLabel',labels18(ord),'XTickLabelRotation',45,...
        'YTick',1:length(ord),'YTickLabel',labels18(ord),'Ydir','normal');
title(sprintf('Population agreement (singular-only), group mean (c=%.2f)', S.coph_r));

% Plot SEM
Ase = S.J_sem(ord,ord);
subplot(1,2,2);
mx = max(Ase(:),[],'omitnan'); if isempty(mx) || mx==0, mx=0.05; end
h2 = imagesc(Ase, [0 mx]); axis square; colormap(parula); colorbar
set(h2,'AlphaData',~isnan(Ase)); set(gca,'Color',[0.9 0.9 0.9],'Ydir','normal');
set(gca,'XTick',[],'YTick',[]); title('SEM across animals');

% (Optional) Dendrogram using same order
hf = figure(2000);clf
% pause(2)
set(hf,'Color','w','Position',[100 200 900 260]);
set(hf,'Units','inches');set(hf,'Position',[1 1 3.5 1.25])
[Hd,~,~] = dendrogram(S.tree, 0, 'Reorder', ord, 'Orientation', 'top','ColorThreshold',0.96); ha = gca;
set(Hd,'LineWidth',0.51); xlim([0.5 26.5]);
set(gca,'XTick',1:length(ord),'XTickLabel',labels18(ord),'XTickLabelRotation',45);
ylabel('Jaccard Distance'); htit = title(sprintf('Hierarchical clustering (c = %.2f)', S.coph_r)); set(htit,'FontWeight','Normal')
changePosition(ha,[-0.0061 0 0.06 0])
format_axes(ha)
save_pdf(hf,mData.pdf_folder,sprintf('OI_Map_cluster.pdf'),600);
%%
A = size(idv,1); B = 2000; cboot = nan(B,1);
for b=1:B
    idx = randi(A,A,1);                      % resample animals
    Jb  = mean(S.J_animal(:,:,idx),3,'omitnan');   % S from pop18_agreement_cluster
    D   = 1 - Jb; D(1:19:end)=0; D(isnan(D))=0;
    Y   = squareform(D,'tovector'); tree_b = linkage(Y,'average');
    cboot(b) = cophenet(tree_b,Y);
end
ci = prctile(cboot,[2.5 97.5]) % 95% CI


    %%
    % === Set metric ===
idv = idvPC;   % later: idv = idvMI;

% Optional nice labels for the 6 cases:
case_labels = {'C3-A1','C3-A2','C4-A1','C4-A2','C5-A1','C5-A2'};

% Build 18 singular-only populations and compute agreement + clustering
[pop, labels18] = build_pop18_singular(idv, case_labels);
labels22 = [{'C2-A1-T','C2-A2-T','C7-A1-T','C7-A2-T'} labels18];

cells_Ab = propsTB.newPC.cells_speed;
pop.X = [cells_Ab cell_popsPC]; 

% cells_Ab = propsTB.newMI.cells_speed;
% pop.X = [cells_Ab cell_popsMI]; 


pop.M = 22;
S = pop18_agreement_cluster(pop);
ord = S.leafOrder;
% ord = 1:18;

% Plot group mean agreement (diag masked), with NaNs as light gray
A18 = S.J_mean(ord,ord);
hf = figure(1000);
set(hf,'Color','w','Position',[100 100 900 360]);
subplot(1,2,1);
h = imagesc(A18, [0 0.3]); axis square; colormap(parula); colorbar
set(h,'AlphaData',~isnan(A18)); set(gca,'Color',[0.9 0.9 0.9]);
set(gca,'XTick',1:22,'XTickLabel',labels22(ord),'XTickLabelRotation',45,...
        'YTick',1:22,'YTickLabel',labels22(ord),'Ydir','normal');
title(sprintf('Population agreement (singular-only), group mean (c=%.2f)', S.coph_r));

% Plot SEM
Ase = S.J_sem(ord,ord);
subplot(1,2,2);
mx = max(Ase(:),[],'omitnan'); if isempty(mx) || mx==0, mx=0.05; end
h2 = imagesc(Ase, [0 mx]); axis square; colormap(parula); colorbar
set(h2,'AlphaData',~isnan(Ase)); set(gca,'Color',[0.9 0.9 0.9],'Ydir','normal');
set(gca,'XTick',[],'YTick',[]); title('SEM across animals');

% (Optional) Dendrogram using same order
hf = figure(2000);clf
% pause(2)
set(hf,'Color','w','Position',[100 200 900 260]);
set(hf,'Units','inches');set(hf,'Position',[1 1 3.5 1.25])
[Hd,~,~] = dendrogram(S.tree, 0, 'Reorder', ord, 'Orientation', 'top','ColorThreshold',0.96); ha = gca;
set(Hd,'LineWidth',0.51); xlim([0.5 18.5]);
set(gca,'XTick',1:22,'XTickLabel',labels22(ord),'XTickLabelRotation',45);
ylabel('Jaccard Distance'); htit = title(sprintf('Hierarchical clustering (c = %.2f)', S.coph_r)); set(htit,'FontWeight','Normal')
changePosition(ha,[-0.0061 0 0.06 0])
format_axes(ha)
save_pdf(hf,mData.pdf_folder,sprintf('OI_Map_cluster.pdf'),600);
