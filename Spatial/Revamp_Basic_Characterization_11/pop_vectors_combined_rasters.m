% combine rasters to show cells' response everywhere. trial wise. in
% reality the stimuli are sequential on top of sequence of trials
%%
si = [Lb_T Ab_T Ar_t_D Ar_i_T ArL_t_D ArL_i_T Ars_t_D Ars_i_T Lbs_T Abs_T ArL_L_T];
Rs = o.Rs(:,si);
ntrials = 50; props1 = get_props_Rs(Rs,ntrials);
good_FR = props1.vals;
good_FR_any = cell_list_op(good_FR,[],'or',1);

[cRs,xss,xse] = combine_rasters(Rs);
mcR = calc_mean_rasters(cRs,[]); 
[CRc,aCRc,mRR] = find_population_vector_corr(cRs,mcR,good_FR_any,0);
%%
an = 1;
ff = makeFigureRowsCols(107,[1 0.5 4 0.5],'RowsCols',[2 1],...
        'spaceRowsCols',[0.05 0.01],'rightUpShifts',[0.11 0.13],'widthHeightAdjustment',...
        [-400 -100]);    set(gcf,'color','w');    set(gcf,'Position',[10 3 6.9 6]);
%     ff = show_population_vector_and_corr(mData,ff,Rs(an,:),mRR(an,:),CRc(an,:),[],[]);
    axes(ff.h_axes(1,1));
    hold on;
    mR = mRR{an,:};
    szy = size(mR,1); szx = size(mR,2);
    imagesc(mR);
    for ii = 1:length(xse)
        plot([xse(ii) xse(ii)],[0 szy],'w');
    end
    xlim([1 szx]); ylim([1 szy]);
    colormap jet
    
    axes(ff.h_axes(2,1)); imagesc(CRc{an,:}); hold on;
    set(gca,'Ydir','normal');
    xlim([1 szx]);
    axis equal
    for ii = 1:length(xse)
        plot([xse(ii) xse(ii)],[0 szx],'w');
        plot([0 szx],[xse(ii) xse(ii)],'w');
    end
    colormap jet
    save_pdf(ff.hf,mData.pdf_folder,sprintf('PV_comb.pdf'),600);