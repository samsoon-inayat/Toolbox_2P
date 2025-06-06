function [mOI,semOI] = heatmap_conj_comp(ha,allresp,type,options)
mData = evalin('base','mData');
siG = options{1};
rasterNamesTxt = options{2};
try
    dis = options{3};
catch
    dis = 0;
end
    
% G = options{3};
[OIo,mOI,semOI,OI_mato,p_vals,h_vals,all_CI,mCI,semCI,all_CI_mat,uni] = get_overlap_index(allresp,0.5,0.05);
    if ~strcmp(get(ha,'type'),'figure')
        axes(ha); hold on;
    else
        ha = gca;
    end
    mUni = nanmean(uni,3); mUni1 = tril(mUni,-1) + tril(mUni,-1)'; mUni2 = triu(mUni,1) + triu(mUni,1)'; mmUni = min(mUni(:)); MmUni = max(mUni(:));
    semUni = nanstd(uni,[],3); semUni1 = tril(semUni,-1) + tril(semUni,-1)'; semUni2 = triu(semUni,1) + triu(semUni,1)'; msemUni = min(semUni(:)); MsemUni = max(semUni(:));
%     mOI = mCI; semOI = semCI;
%     mSel = mCI;
%     mOI = mSel;
    switch type
        case 0
            n = 0;
            mOI = mOI; semOI = semOI;
        case 1
            mOI = mCI; semOI = semCI;
        case 2
            mOI = mUni1; semOI = semUni1;
        case 3
            mOI = mUni2; semOI = semUni2;
    end
    sz = size(mOI,1);
%     mOI = OI_mat(:,:,4);
    oM = ones(size(mOI));
    mask1 = (triu(oM,0) & tril(oM,0)); mOI(mask1==1) = NaN;
    maxI = max([mOI(:);semOI(:)]);
    minI = min([mOI(:);semOI(:)]);
%     minI = 0; maxI = 0.6;
    mask = tril(NaN(size(mOI)),0); mask(mask==0) = 1; 
%     mOI = mOI .* mask;
    imAlpha=ones(size(mOI));    %imAlpha(isnan(mask))=0.25; 
    imAlpha(mask1 == 1) = 0;
%     ff = makeFigureRowsCols(2020,[10 4 6 1.5],'RowsCols',[1 2],'spaceRowsCols',[0.1 0.01],'rightUpShifts',[0.05 0.8],'widthHeightAdjustment',[-240 -150]);
%     hf = get_figure(6,[8 3 3.5 3.5]);
%     hf = get_figure(6,[8 3 4 4]);
    %
%     axes(ff.h_axes(1));
    im1 = imagesc(mOI,[minI,maxI]);    im1.AlphaData = imAlpha;
    if ~dis
    for ii = 1:18
        plot([(10.5+((ii-1)*10)) (10.5+((ii-1)*10))],[0 180.5],'w','linewidth',0.1); 
        plot([0 180.5],[(10.5+((ii-1)*10)) (10.5+((ii-1)*10))],'w','linewidth',0.1); 
    end
    end
%     for ii = [2 6 9 12]   
%         plot([(10.5+((ii-1)*10)) (10.5+((ii-1)*10))],[0 80.5],'k','linewidth',1); 
%         plot([0 80.5],[(10.5+((ii-1)*10)) (10.5+((ii-1)*10))],'k','linewidth',1); 
%     end
    set(gca,'color',0.5*[1 1 1]);    colormap parula;    %axis equal
    format_axes(gca);
    set_axes_limits(gca,[0.5 sz+0.5],[0.5 sz+0.5]);
    ttxl = rasterNamesTxt(siG);
    if dis
        xtickvals = 1:size(mOI,2);%[5 15 25 60 100 115 125];
    else
        xtickvals = 5:10:size(mOI,2);%[5 15 25 60 100 115 125];
    end
    xticklabels = rasterNamesTxt(siG);

    set(gca,'xtick',xtickvals,'ytick',xtickvals,'xticklabels',xticklabels,'yticklabels',xticklabels,'Ydir','normal'); xtickangle(45);%ytickangle(45);
    yyaxis right
    set(gca,'ytick',xtickvals,'yticklabels',yticklabels,'tickdir','out');
    box off
    changePosition(gca,[-0.01 0.00 0.037 0.033]);
    hc = putColorBar(gca,[0.1 -0.08 -0.2 0.03],{sprintf('%.1f',minI),sprintf('%.1f',maxI)},6,'northoutside',[0.07 0.09 0.02 0.09]);
    colormap jet
%     save_pdf(hf,mData.pdf_folder,sprintf('OI_Map_mean_tu_spatial_%c.pdf',G),600);