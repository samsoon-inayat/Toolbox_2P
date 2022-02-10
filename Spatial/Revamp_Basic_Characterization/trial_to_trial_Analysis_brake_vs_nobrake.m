function trial_to_trial_Analysis_brake_vs_nobrake

%% find spatial trial to trial correlation
while 1
    trialNums = [1:10];
   si = [Lb Ab_On Ab_Off ArL_L Ar_On ArL_On Ars_On Ar_Off ArL_Off Ars_Off Lbs Abs_On Abs_Off ];
    Rs = o.Rs(:,si);mR = o.mR(:,si); RsG = Rs; siG = si; propsG = get_props_Rs(RsG,[40,100]); respG = propsG.vals;
    avgProps = get_props_Rs(Rs,[40,100]); respM = avgProps.good_FR;
    for cn = 1:length(si)
        trials = mat2cell([1:10]',ones(size([1:10]')));
        trials = mat2cell([trialNums]',ones(size([trialNums]')));
        RsC = repmat(Rs(:,cn),1,10);
        mRsCT = cell(size(RsC,1),length(trials));
        for ii = 1:length(trials)
            ii;
            [mRsCT(:,ii),~] = calc_mean_rasters(RsC(:,1),trials{ii});
        end
        allmRsT{cn} = mRsCT;
        allRsC{cn} = RsC;
    end
    disp('Done');
    %%
    break;
end

%% find remapping (code will take long time to run)
while 1
    si = [Ar_t_D ArL_t_D Ars_t_D];
    si = [Ar_t_D Ar_i_T ArL_t_D ArL_i_T Ars_t_D Ars_i_T];
    si = [Lb_T ArL_L_T Lbs_T Ab_T Abs_T Ar_t_D ArL_t_D Ars_t_D Ar_i_T ArL_i_T Ars_i_T];
%     si = [Lb_T Ab_T Ar_t_D Ar_i_T ArL_t_D ArL_i_T Ars_t_D Ars_i_T Lbs_T Abs_T];
%     si = [Ar_i_T ArL_i_T Ars_i_T];
    Rs = o.Rs(:,si);mR = o.mR(:,si); RsG = Rs;
    avgProps = get_props_Rs(Rs,[0,40]); respG{1} = avgProps.good_FR; avgProps = get_props_Rs(Rs,[40,70]); respG{2} = avgProps.good_FR;
    avgProps = get_props_Rs(Rs,[70,100]); respG{3} = avgProps.good_FR; avgProps = get_props_Rs(Rs,[40,100]); respG{4} = avgProps.good_FR;
    for cn = 7%:length(si)
        trials = mat2cell([1:10]',ones(size([1:10]')));
        RsC = repmat(Rs(:,cn),1,10);
        RsC = find_responsive_rasters(RsC,1:10);
        mRsCT = allmRsT{cn};
        for pn = 1:4
            [cn pn]
            out_C{cn,pn} = find_population_vector_corr_remap(RsC,mRsCT,respG{pn});
        end
    end
    %%
    break;
end

%% average correlation of all animals
while 1
    trialsRemap = out_C{7,1};
    meanPVcorr = trialsRemap.mean_PV_corr;
    bigMeanPVCorr = [];
    for rr = 1:size(meanPVcorr,1)
        tempBM = [];
        for cc = 1:size(meanPVcorr,2)
            tempBM = [tempBM meanPVcorr{rr,cc}];
        end
        bigMeanPVCorr = [bigMeanPVCorr;tempBM];
    end
    hf = get_figure(6,[8 3 7 7]);
    imagesc(bigMeanPVCorr);
    %
    %     axes(ff.h_axes(1));
    plot([10.5 10.5],[0 30.5],'r'); plot([20.5 20.5],[0 30.5],'r');
    plot([0 30.5],[10.5 10.5],'r'); plot([0 30.5],[20.5 20.5],'r');
    set(gca,'color',0.5*[1 1 1]);    colormap parula;    %axis equal
    format_axes(gca);
    set_axes_limits(gca,[0.5 sz+0.5],[0.5 sz+0.5]);
    set(gca,'xtick',1:length(txl),'ytick',1:length(txl),'xticklabels',txl,'yticklabels',txl,'Ydir','reverse'); xtickangle(45);
    text(-0.3,sz+1.1,'Average Overlap Index (N = 5 animals)','FontSize',5); set(gca,'Ydir','normal');ytickangle(20);
    box on
    changePosition(gca,[0.03 0 -0.04 0]);
    hc = putColorBar(gca,[0.09 0.07 -0.11 -0.15],{sprintf('%d',minI),sprintf('%.1f',maxI)},6,'eastoutside',[0.1 0.05 0.06 0.05]);
    
%     ff = makeFigureRowsCols(106,[1 0.5 9 9],'RowsCols',[10 10],...
%         'spaceRowsCols',[0.02 0.02],'rightUpShifts',[0.1 0.1],'widthHeightAdjustment',...
%         [-40 -40]);
%     set(gcf,'color','w');
%     set(gcf,'Position',[5 1 9 9]);
%     ff = show_remapping_corr_plots(mData,ff,trialsRemap.mean_PV_corr,trialsRemap.xs,[]);
%     save_pdf(ff.hf,mData.pdf_folder,sprintf('remap_corr_C.pdf'),600);
    %%
    break;
end

%% Overlap Indices ImageSC all
while 1
    avgProps = get_props_Rs(RsG,[50,100]); 
    respG = avgProps.good_FR_and_untuned;
    an  = 1:5; eic = 1; sp = 0; intersect_with_global = 0; only_global = 0;
    allresp = []; ind = 1;
    all_peakL = [];
    for cn = 1:length(si)
        mRsCT = allmRsT{cn};
        resp = []; peak_locations = [];
        for rr = 1:size(mRsCT,1)
            for cc = 1:size(mRsCT,2)
                this_mat = mRsCT{rr,cc};
                [~,peakL] = max(this_mat,[],2);
%                 size_tmat(rr,cc) = size(this_mat,2);
                resp{rr,cc} = sum(this_mat,2) > 0;
                if intersect_with_global
                    resp{rr,cc} = resp{rr,cc} & respG{rr,cn};
                end
                if only_global
                    resp{rr,cc} = respG{rr,cn};
                end
                if sp == 1
                    if cn == 1
                        respSe = respSeL{eic}; resp{rr,cc} = resp{rr,cc} & respSe{rr,1};
                    end
                    if cn == 2
                        respSe = respSeL{eic}; resp{rr,cc} = resp{rr,cc} & respSe{rr,2};
                    end
                    if cn == 3
                        respSe = respSeL{eic}; resp{rr,cc} = resp{rr,cc} & respSe{rr,3};
                    end
                    if cn == 4
                        respSe = respSeA{eic}; resp{rr,cc} = resp{rr,cc} & respSe{rr,1};
                    end
                    if cn == 5
                        respSe = respSeA{eic}; resp{rr,cc} = resp{rr,cc} & respSe{rr,2};
                    end
                end
                peakL(~resp{rr,cc}) = NaN;
                peak_locations{rr,cc} = peakL;
                if rr == 1
                    txl{ind} = sprintf('C%dT%d',cn,cc);
                    ind = ind + 1;
                end
            end
%             oc(rr,cn) = find_cells_based_on_cluster(cell2mat(resp(rr,:)));
        end
        allresp = [allresp resp]; all_peakL = [all_peakL peak_locations];
        
    end
    i_allresp = cell_list_op(allresp,[],'not');
    [OI,mOI,semOI,OI_mat,p_vals,h_vals] = get_overlap_index(allresp(an,:),0.5,0.05,0);
%     break;
% end
% %%
% 
% %%
% while 1   
%     [OI,mOI,semOI,OI_mat,p_vals,h_vals] = get_overlap_index(allresp,0.5,0.05,0);
%     mOI = OI{4}; semOI = semOI;
    
%     [OI,mOI,semOI,OI_mat,p_vals,h_vals] = get_overlap_index(i_allresp(an,:),0.5,0.05);
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
%     ff = makeFigureRowsCols(2020,[10 4 6 1.5],'RowsCols',[1 2],'spaceRowsCols',[0.1 0.01],'rightUpShifts',[0.05 0.13],'widthHeightAdjustment',[-240 -150]);
    hf = get_figure(6,[8 3 3.5 3.5]);
%     hf = get_figure(6,[8 3 4 4]);
    %
%     axes(ff.h_axes(1));
    im1 = imagesc(mOI,[minI,maxI]);    im1.AlphaData = imAlpha;
    
    for ii = 1:12
        plot([(10.5+((ii-1)*10)) (10.5+((ii-1)*10))],[0 130.5],'w','linewidth',0.1); 
        plot([0 130.5],[(10.5+((ii-1)*10)) (10.5+((ii-1)*10))],'w','linewidth',0.1); 
    end
    for ii = [3 10]   
        plot([(10.5+((ii-1)*10)) (10.5+((ii-1)*10))],[0 130.5],'k','linewidth',1); 
        plot([0 130.5],[(10.5+((ii-1)*10)) (10.5+((ii-1)*10))],'k','linewidth',1); 
    end
    set(gca,'color',0.5*[1 1 1]);    colormap parula;    %axis equal
    format_axes(gca);
    set_axes_limits(gca,[0.5 sz+0.5],[0.5 sz+0.5]);
    ttxl = rasterNamesTxt(siG);
    xtickvals = 5:10:130;%[5 15 25 60 100 115 125];
    % si = [Lb Ab_On Ab_Off ArL_L Ar_On ArL_On Ars_On Ar_Off ArL_Off Ars_Off Lbs Abs_On Abs_Off ];
    xticklabels = {'Lb - ON','Lb* - ON','Ab - ON','Ab* - ON','Ab - OFF','Ab* - OFF','ArL-L - ON','Ar - ON','ArL - ON','Ar* - ON','Ar - OFF','ArL - OFF','Ar* - OFF',};
    xticklabels = {'Lb - ON','Ab - ON','Ab - OFF','ArL-L - ON','Ar - ON','ArL - ON','Ar* - ON','Ar - OFF','ArL - OFF','Ar* - OFF','Lb* - ON','Ab* - ON','Ab* - OFF'};
    yticklabels = {'ON','ON','ON','ON','ON','OFF','OFF','ON','ON','ON','OFF','OFF','OFF',};

    set(gca,'xtick',xtickvals,'ytick',xtickvals,'xticklabels',xticklabels,'yticklabels',xticklabels,'Ydir','normal'); xtickangle(45);%ytickangle(45);
    yyaxis right
    set(gca,'ytick',xtickvals,'yticklabels',yticklabels,'tickdir','out');
%     text(-0.3,sz+1.1,'Average Overlap Index (N = 5 animals)','FontSize',5); set(gca,'Ydir','normal');ytickangle(20);
    box off
    changePosition(gca,[0.03 0.007 0.037 0.033]);
    hc = putColorBar(gca,[0.15 -0.08 -0.2 0.03],{sprintf('%d',minI),sprintf('%.1f',maxI)},6,'northoutside',[0.03 0.09 0.02 0.09]);
    colormap jet
    save_pdf(hf,mData.pdf_folder,sprintf('OI_Map_mean_tu_spatial.pdf'),600);
    %%
    break;
end


