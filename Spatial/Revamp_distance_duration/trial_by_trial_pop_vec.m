function trial_by_trial_pop_vec





%% population vector and correlation sensory
while 1
    an = 1; cnseq = 6; 
    tRs = [RsG(an,cn) allRsC{cn}(an,:)]; tmR = [mRsG(an,cn) allmRsT{cn}(an,:)];
    
%     respT = FD_conj{1}(an,1);
%     respT = FT_Dur_comp{1}(an,1);
    cn = 1;
%     respT = cell_list_op(respG(an,cn),respG(an,cn+3),'or');
    respT = respG(an,cn);
    respT = allresp_OR(an,cn);
%     respC = dis_cells_T;%cell_list_op(dis_cells_T,dur_cells_I,'and');
%     respC = dur_cells_I;
%     respT = respC(an,3);
%     respT = allresp_OR;
    [~,~,seq] = findPopulationVectorPlot(mRsG{an,1},respT{1});
%     [~,~,seq] = findPopulationVectorPlot(allmRsT{cn}{an,1},respT{1});
    
    cn = 1;
    tRs = [RsG(an,cn) allRsC{cn}(an,:)]; tmR = [allmRsT{cn}(an,:) mRsG(an,cn)];
    [ff] = plot_pop_vec_trials(100,[0.5 3 19 5.5],tRs,tmR,respT,seq);
    
%     tRs = [RsG(an,cn+3) allRsC{cn+3}(an,:)]; tmR = [allmRsT{cn+3}(an,:) mRsG(an,cn+3)];
%     ff = plot_pop_vec_trials(101,[10 3 9 4.5],tRs,tmR,respT,seq)
    
%     cn = 2;
%     tRs = [RsG(an,cn) allRsC{cn}(an,:)]; tmR = [mRsG(an,cn) allmRsT{cn}(an,:)];
%     [ff] = plot_pop_vec_trials(102,[0.5 5 7 1.5],tRs,tmR,respT,seq)
%     
%     tRs = [RsG(an,cn+3) allRsC{cn+3}(an,:)]; tmR = [mRsG(an,cn+3) allmRsT{cn+3}(an,:)];
%     ff = plot_pop_vec_trials(103,[8 5 7 1.5],tRs,tmR,respT,seq)
%     
%     cn = 3;
%     tRs = [RsG(an,cn) allRsC{cn}(an,:)]; tmR = [mRsG(an,cn) allmRsT{cn}(an,:)];
%     [ff] = plot_pop_vec_trials(104,[0.5 7 7 1.5],tRs,tmR,respT,seq)
%     
%     tRs = [RsG(an,cn+3) allRsC{cn+3}(an,:)]; tmR = [mRsG(an,cn+3) allmRsT{cn+3}(an,:)];
%     ff = plot_pop_vec_trials(105,[8 7 7 1.5],tRs,tmR,respT,seq)
    

    %%
    cn = cns(2);
% %     si = [Lb Ab_On Ab_Off ArL_L Ar_On ArL_On Ars_On Ar_Off ArL_Off Ars_Off Lbs];
%     props1 = get_props_Rs(o.Rs(:,si),50);
    resp = gFR_T_g_D;
%     resp = props1.good_FR_and_untuned;
%     Rs = [o.Rs(an,si(cn)) allRsC{cn}(an,:)];mR = [o.mR(an,si(cn)) allmRsT{cn}(an,:)];
    respT = repmat(resp(an,cn),1,size(mR,2));
    for cc = 1:size(mR,2)
        this_mat = mR{1,cc};
    end
    ff = makeFigureRowsCols(106,[1 0.5 4 0.5],'RowsCols',[2 11],...
        'spaceRowsCols',[0.02 0.025],'rightUpShifts',[0.03 0.13],'widthHeightAdjustment',...
        [-30 -200]);    set(gcf,'color','w');    set(gcf,'Position',[1 4 14 3]);
    [CRc,aCRc,mRR] = find_population_vector_corr(Rs,mR,respT,1);
    ff = show_population_vector_and_corr(mData,ff,Rs,mRR,CRc,[],[]);
    changePosition(ff.h_axes(2,1).YLabel,[-2.5 0 0]); 
    colormap parula
%     save_pdf(ff.hf,mData.pdf_folder,sprintf('PV_trials_%s.pdf',selected_property),600);
    %%
    an = 1; cn = 5;%cns(3);
    tn = 5;
%     si = [Lb Ab_On Ab_Off ArL_L Ar_On ArL_On Ars_On Ar_Off ArL_Off Ars_Off Lbs];
%     props1 = get_props_Rs(o.Rs(:,si),50);
%     resp = props1.inh;
    resp = props1.good_FR;
    Rs = [RsG(an,cn) allRsC{cn}(an,:)];mR = [mRsG(an,cn) allmRsT{cn}(an,:)];%mR = [o.mR(an,si(cn)) allmRsT{cn}(an,:)];
    respPV = allresp_trials{an,cn}; for tni = 1:10 respPVC(tni) = {respPV(:,tni)}; end
    respT = repmat({respPV(:,tn)},1,size(mR,2));%respT = repmat(respG(an,cn),1,size(mR,2));
    respT = [respG(an,cn),respPVC];
    for cc = 1:size(mR,2)
        this_mat = mR{1,cc};
    end
    ff = makeFigureRowsCols(107,[1 0.5 4 0.5],'RowsCols',[2 11],...
        'spaceRowsCols',[0.02 0.025],'rightUpShifts',[0.03 0.13],'widthHeightAdjustment',...
        [-30 -200]);    set(gcf,'color','w');    set(gcf,'Position',[1 7 14 3]);
%     [CRc,aCRc,mRR] = find_population_vector_corr(Rs,mR,respT,tn+1);
    [CRc,aCRc,mRR] = find_population_vector_corr(Rs,mR,respT,[]);
    ff = show_population_vector_and_corr(mData,ff,Rs,mRR,CRc,[],[]);
    changePosition(ff.h_axes(2,1).YLabel,[-2.5 0 0]); 
    colormap parula
%     save_pdf(ff.hf,mData.pdf_folder,sprintf('PV_trials_%s.pdf',selected_property),600);
    
    
    
%     save_pdf(ff.hf,mData.pdf_folder,sprintf('PV_trials_%s.pdf',selected_property),600);

    %%
    break;
end
