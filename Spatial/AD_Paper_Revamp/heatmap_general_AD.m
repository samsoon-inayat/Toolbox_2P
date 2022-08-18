
ntrials = 50;
si = [C1_t_D C2_t_D C3_t_D C4_t_D C1_i_T C2_i_T C3_i_T C4_i_T];
%     si = [C1_i_T C2_i_T C3_i_T C4_i_T];
Rs_C = oC.Rs(:,si); Rs_A = oA.Rs(:,si); mRs_C = oC.mR(:,si); mRs_A = oA.mR(:,si);
props_C = get_props_Rs(oC.Rs(:,si),ntrials); props_A = get_props_Rs(oA.Rs(:,si),ntrials);
% pop_var_name = {'all','vals','valsT','Nvals','good_zMI','Ngood_zMI'};
pop_var_name = {'good_zMI','good_Gauss','good_MFR'};
pop_var_name = {'all'};
pop_var_name = {'vals'};
sel_pop_C = cell_list_op(props_C,pop_var_name); sel_pop_A = cell_list_op(props_A,pop_var_name);

%% for heat maps
figdim = 3;
% hf = get_figure(5,[9 2 figdim figdim]);

all_cells_list = sel_pop_C; event_type = {'C1-D','C2-D','C3-D','C4-D','C1-T','C2-T','C3-T','C4-T'};
%     all_cells_list = all_gV_A;
sh = 0;
good_FR = circshift(all_cells_list,sh,2);
txl = circshift(event_type,sh,2);
seq_nums = circshift(1:size(all_cells_list,2),sh,2);

[OIo,mOIo,semOIo,OI_mat,p_vals,h_vals,all_CI,mCI,semCI,all_CI_mat,uni] = get_overlap_index(good_FR,0.5,0.05);
mOI = mCI; semOI = semCI;


mUni = nanmean(uni,3); mUni1 = tril(mUni,-1) + tril(mUni,-1)'; mUni2 = triu(mUni,1) + triu(mUni,1)'; mmUni = min(mUni(:)); MmUni = max(mUni(:));
semUni = nanstd(uni,[],3)/sqrt(5); semUni1 = tril(semUni,-1) + tril(semUni,-1)'; semUni2 = triu(semUni,1) + triu(semUni,1)'; msemUni = min(semUni(:)); MsemUni = max(semUni(:));
%     figure(1000);clf;subplot 131;imagesc(mUni,[mmUni MmUni]);set(gca,'YDir','normal');subplot 132;imagesc(mUni1,[mmUni MmUni]);set(gca,'YDir','normal');subplot 133;imagesc(mUni2,[mmUni MmUni]);set(gca,'YDir','normal');


ff = makeFigureRowsCols(107,[10 3 3.5 1.4],'RowsCols',[1 2],'spaceRowsCols',[0.2 0.2],'rightUpShifts',[0.2 0.17],'widthHeightAdjustment',[-10 -300]);
set(gcf,'color','w');     ylims = [0 1];
stp = 0.22; widths = 0.5*([2 2 2 2 0.4 0.4]-0.2); gap = 0.5*1.75; adjust_axes(ff,ylims,stp,widths,gap,{'Euclidean Distance'});
gap1 = 0.24; widths1 = 0.5*0.75; height1 = widths1 - 0.05;

mats = {mCI,mUni}; semmats = {semCI,semUni};
titles = {'Conjunction','Complementation'};
for ii = 1:length(mats)
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
                text(rr-0.25,cc,sprintf('%.0f',mCI(rr,cc)),'FontSize',5);
            end
        end
    end
%     plot([(10.5+((ii-1)*10)) (10.5+((ii-1)*10))],[0 150.5],'w','linewidth',0.1); 
%     plot([0 150.5],[(10.5+((ii-1)*10)) (10.5+((ii-1)*10))],'w','linewidth',0.1); 
    set(gca,'Ydir','normal');ytickangle(15);        box on
    format_axes_b(gca);
    [hc,hca] = putColorBar(gca,[0.09 0.07 -0.11 -0.15],{sprintf('%.1f',minI),sprintf('%.1f',maxI)},6,'eastoutside',[0.07 0.05 0.07 0.05]);
    colormap jet
    htit = title('Mean'); changePosition(htit,[0 0.07 0]); set(htit,'FontSize',6);
    ht = set_axes_top_text_no_line(gcf,gca,titles{ii},[0 0.05 0 0]);set(ht,'FontSize',6,'FontWeight','Bold')

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
    htit = title('SEM'); changePosition(htit,[0 -0.2 0]); set(htit,'FontSize',6);
    box on
    format_axes_b(gca);
    colormap jet
end
save_pdf(ff.hf,mData.pdf_folder,sprintf('OI_Map_%d_%d.pdf',ntrials,sh),600);
%%
ff = makeFigureRowsCols(107,[10 3 3.5 1.5],'RowsCols',[2 3],'spaceRowsCols',[0.3 -0.02],'rightUpShifts',[1 0.13],'widthHeightAdjustment',[-10 -300]);
set(gcf,'color','w');     ylims = [0 1];
stp = 0.25; widths = [2 2 2 0.4 0.4 0.4]+0.07; gap = 0.13; adjust_axes(ff,ylims,stp,widths,gap,{'Euclidean Distance'});
stp = 0.35; widths = 0.45*([2 2 2 0.4 0.4 0.4]+0.07); gap = 0.15; adjust_axes(ff,ylims,stp,widths,gap,{'Euclidean Distance'});
%
hms = {mCI,mUni1,mUni2};
ahc_col_th = 0.7;

hms1 = {100-mCI,mUni1,mUni2};
for hi = 1:length(hms1)
    mOI1 = hms1{hi}; mOI1(mask1==1) = 0; Di = squareform(mOI1,'tovector');%pdist(mOI1,@naneucdist);
    tree = linkage(mOI1,'average'); [c,d] = cophenet(tree,Di); r = corr(Di',d','type','spearman');    leafOrder = optimalleaforder(tree,Di);
    hf = figure(100000000); %     leafOrder1 = leafOrder([1:3 10:12 4:9]);
%     leafOrder1 = circshift(leafOrder,3);
    figure(hf);clf
    [H,T,TC] = dendrogram(tree,'Orientation','top','ColorThreshold',ahc_col_th*max(tree(:,3)),'Reorder',leafOrder);ylims = ylim; xlims = xlim;    set(gca,'xtick',[],'ytick',[]);
    set(gcf,'units','inches'); set(gcf,'Position',[5 2 0.9 0.5]); close(hf);
    axes(ff.h_axes(1,hi));
    [H,T,TC] = dendrogram(tree,'Orientation','top','ColorThreshold',ahc_col_th*max(tree(:,3)),'Reorder',leafOrder);
    set(H,'linewidth',0.5); set(gca,'xticklabels',txl(leafOrder));xtickangle(45); format_axes(gca);
    if hi == 1
        hx = ylabel({'Euc. Dist.'});changePosition(hx,[0 0 0]);
    end
%     xlim([xlims(1)+0.5 xlims(2)-0.5]);
%         changePosition(gca,[0.0 0.0 0.07 0.05]); text(0.5,ylims(2)+0,sprintf('CC = %.2f',c),'FontSize',6);
    hcct = set_axes_top_text_no_line(ff.hf,gca,sprintf('CC = %.2f',c),[0 -0.05 0 0]); set(hcct,'FontWeight','Normal')
end

for hi = 1:length(hms)
    mOI1 = hms{hi}; mOI1(mask1==1) = NaN; Di = pdist(mOI1,@naneucdist);
    tree = linkage(Di,'average'); [c,d] = cophenet(tree,Di); r = corr(Di',d','type','spearman');    leafOrder = optimalleaforder(tree,Di);
    hf = figure(100000000); %     leafOrder1 = leafOrder([1:3 10:12 4:9]);
%     leafOrder1 = circshift(leafOrder,3);
    figure(hf);clf
    [H,T,TC] = dendrogram(tree,'Orientation','top','ColorThreshold',ahc_col_th*max(tree(:,3)),'Reorder',leafOrder);ylims = ylim; xlims = xlim;    set(gca,'xtick',[],'ytick',[]);
    set(gcf,'units','inches'); set(gcf,'Position',[5 2 0.9 0.5]); close(hf);
    axes(ff.h_axes(2,hi));
    [H,T,TC] = dendrogram(tree,'Orientation','top','ColorThreshold',ahc_col_th*max(tree(:,3)),'Reorder',leafOrder);
    set(H,'linewidth',0.5); set(gca,'xticklabels',txl(leafOrder));xtickangle(45); format_axes(gca);
    if hi == 1
        hx = ylabel({'Euc. Dist.'});changePosition(hx,[0 0 0]);
    end
%     xlim([xlims(1)+0.5 xlims(2)-0.5]);
%         changePosition(gca,[0.0 0.0 0.07 0.05]); text(0.5,ylims(2)+0,sprintf('CC = %.2f',c),'FontSize',6);
    hcct = set_axes_top_text_no_line(ff.hf,gca,sprintf('CC = %.2f',c),[0 0.05 0 0]); set(hcct,'FontWeight','Normal')
%     set_axes_top_text(ff.hf,ff.h_axes(1),sprintf('Cophenetic Correlation = %.2f',c));
end
save_pdf(ff.hf,mData.pdf_folder,sprintf('OI_Map_cluster_%d.pdf',sh),600);

%%
ff = makeFigureRowsCols(107,[10 3 3.5 1],'RowsCols',[1 3],'spaceRowsCols',[0.3 -0.02],'rightUpShifts',[1 0.23],'widthHeightAdjustment',[-10 -550]);
set(gcf,'color','w');     ylims = [0 1];
% stp = 0.25; widths = [2 2 2 0.4 0.4 0.4]+0.07; gap = 0.13; adjust_axes(ff,ylims,stp,widths,gap,{'Euclidean Distance'});
stp = 0.25; widths = 0.45*([2 2 2 0.4 0.4 0.4]+0.07); gap = 0.15; adjust_axes(ff,ylims,stp,widths,gap,{'Euclidean Distance'});
%
hms = {mCI,mUni1,mUni2};
ahc_col_th = 0.7;


for hi = 1:length(hms)
    mOI1 = hms{hi}; mOI1(mask1==1) = NaN; Di = pdist(mOI1,@naneucdist);
    tree = linkage(Di,'average'); [c,d] = cophenet(tree,Di); r = corr(Di',d','type','spearman');    leafOrder = optimalleaforder(tree,Di);
    hf = figure(100000000); %     leafOrder1 = leafOrder([1:3 10:12 4:9]);
%     leafOrder1 = circshift(leafOrder,3);
    figure(hf);clf
    [H,T,TC] = dendrogram(tree,'Orientation','top','ColorThreshold',ahc_col_th*max(tree(:,3)),'Reorder',leafOrder);ylims = ylim; xlims = xlim;    set(gca,'xtick',[],'ytick',[]);
    set(gcf,'units','inches'); set(gcf,'Position',[5 2 0.9 0.5]); close(hf);
    axes(ff.h_axes(1,hi));
    [H,T,TC] = dendrogram(tree,'Orientation','top','ColorThreshold',ahc_col_th*max(tree(:,3)),'Reorder',leafOrder);
    set(H,'linewidth',0.5); set(gca,'xticklabels',txl(leafOrder));xtickangle(45); format_axes(gca);
    if hi == 1
        hx = ylabel({'Euc. Dist.'});changePosition(hx,[0 0 0]);
    end
%     xlim([xlims(1)+0.5 xlims(2)-0.5]);
%         changePosition(gca,[0.0 0.0 0.07 0.05]); text(0.5,ylims(2)+0,sprintf('CC = %.2f',c),'FontSize',6);
    hcct = set_axes_top_text_no_line(ff.hf,gca,sprintf('CC = %.2f',c),[0 0.05 0 0]); set(hcct,'FontWeight','Normal')
%     set_axes_top_text(ff.hf,ff.h_axes(1),sprintf('Cophenetic Correlation = %.2f',c));
end
save_pdf(ff.hf,mData.pdf_folder,sprintf('OI_Map_cluster_%d.pdf',sh),600);


