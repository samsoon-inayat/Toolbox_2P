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
    an = 4; cn = 9; tn = 1;
    resp = props1.good_FR;
    Rs = [RsG(an,cn) allRsC{cn}(an,:)];mR = [mRsG(an,cn) allmRsT{cn}(an,:)];%mR = [o.mR(an,si(cn)) allmRsT{cn}(an,:)];
    respPV = allresp_trials{an,cn}; for tni = 1:10 respPVC(tni) = {respPV(:,tni)}; end
    respT = [respG(an,cn),respPVC];
%     respT = repmat(respG(an,cn),1,size(mR,2));
    for cc = 1:size(mR,2)
        this_mat = mR{1,cc};
    end
    ff = makeFigureRowsCols(108,[1.5 7 14 2.25],'RowsCols',[2 11],...
        'spaceRowsCols',[0.05 0.025],'rightUpShifts',[0.03 0.13],'widthHeightAdjustment',...
        [-30 -120]);   
    [CRc,aCRc,mRR] = find_population_vector_corr(Rs,mR,respT,[]);
    ff = show_population_vector_and_corr(mData,ff,Rs,mRR,CRc,[],[]);
    changePosition(ff.h_axes(2,1).YLabel,[-2.5 0 0]); 
    colormap parula

    %%
    break;
end



%%
    an = 4; cn = 3; tn = 1;
    resp = props1.good_FR;
    Rs = [RsG(an,cn) allRsC{cn}(an,:)];mR = [mRsG(an,cn) allmRsT{cn}(an,:)];%mR = [o.mR(an,si(cn)) allmRsT{cn}(an,:)];
    respPV = allresp_trials{an,cn}; for tni = 1:10 respPVC(tni) = {respPV(:,tni)}; end
    respT = [respG(an,cn),respPVC];
%     respT = repmat(respG(an,cn),1,size(mR,2));
    for cc = 1:size(mR,2)
        this_mat = mR{1,cc};
    end
    ff = makeFigureRowsCols(108,[1.5 7 14 2.25],'RowsCols',[2 11],...
        'spaceRowsCols',[0.05 0.025],'rightUpShifts',[0.03 0.13],'widthHeightAdjustment',...
        [-30 -120]);   
    [CRc,aCRc,mRR] = find_population_vector_corr(Rs,mR,respT,[]);
    ff = show_population_vector_and_corr(mData,ff,Rs,mRR,CRc,[],[]);
    changePosition(ff.h_axes(2,1).YLabel,[-2.5 0 0]); 
    colormap parula
    
    
    
%     an = 4; cn = 6;%cns(3);
%     tn = 1;
    resp = props1.good_FR;
    Rs = [RsG(an,cn) allRsC{cn}(an,:)];mR = [mRsG(an,cn) allmRsT{cn}(an,:)];%mR = [o.mR(an,si(cn)) allmRsT{cn}(an,:)];
    respPV = allresp_trials{an,cn}; for tni = 1:10 respPVC(tni) = {respPV(:,tni)}; end
    respT = repmat({respPV(:,tn) & respPV(:,tn+1)},1,size(mR,2));%respT = repmat(respG(an,cn),1,size(mR,2));

    for cc = 1:size(mR,2)
        this_mat = mR{1,cc};
    end
      ff = makeFigureRowsCols(107,[1.5 4 14 2.25],'RowsCols',[2 11],...
        'spaceRowsCols',[0.05 0.025],'rightUpShifts',[0.03 0.13],'widthHeightAdjustment',...
        [-30 -120]);    
    [CRc,aCRc,mRR] = find_population_vector_corr(Rs,mR,respT,tn+1);
%     [CRc,aCRc,mRR] = find_population_vector_corr(Rs,mR,respT,[]);
    ff = show_population_vector_and_corr(mData,ff,Rs,mRR,CRc,[],[]);
    changePosition(ff.h_axes(2,1).YLabel,[-2.5 0 0]); 
    colormap parula
    %%
    cns = [1 2 9 10];
    ff = makeFigureRowsCols(108,[1 1 16 9],'RowsCols',[2 length(cns)],...
    'spaceRowsCols',[0.05 0.025],'rightUpShifts',[0.03 0.13],'widthHeightAdjustment',...
    [-30 -120]);   

    an = 1; 
    
    for cni = 1:length(cns)
        cn = cns(cni);
        Rs = [RsG(an,cn) allRsC{cn}(an,:)];mR = [mRsG(an,cn) allmRsT{cn}(an,:)];%mR = [o.mR(an,si(cn)) allmRsT{cn}(an,:)];
        respPV = allresp_trialsC{an,cn};    respPVor = cell_list_op(respPV,[],'or',1);
        mR = allmRsT{cn}(an,:);
        respT = respG(an,cn);
        respT = respPV;
        respT = respPVor;
%         respT = resp_All(an,cn);
        [pv,pvc,pvct,pvg] = get_population_vector(mR,respT);
        axes(ff.h_axes(1,cni));imagesc(pvg.pv); colorbar; set(gca,'YDir','Normal');
        axes(ff.h_axes(2,cni));imagesc(pvg.pvc); colorbar; set(gca,'YDir','Normal');
    end
%     hf = get_figure(8,[5 1 9 9]);hold on;
%     imagesc(pvg.pvc); colorbar
    %%
    otd = []; all_otd = [];
    for an = 1:5
        tnii = 1;
        for cn = 1:length(si)
            respPV = allresp_trialsC{an,cn};
            for tni = 2:10
               otd(tni-1,cn,an) = sum(sum([respPV{tni-1} & ~respPV{tni} ~respPV{tni-1} & respPV{tni}]))/size(respPV{1},1);
               all_otd(an,tnii) = otd(tni-1,cn,an); tnii = tnii + 1;
            end
        end
    end
    m_otd_tr = squeeze(mean(otd,1))';
   
    inds = [1 2 3 4];
    varC = all_otd(:,1:36);
    [within,dvn,xlabels,awithinD] = make_within_table({'Co','Ph','Tr'},[2,2,9]);
    dataT = make_between_table({varC},dvn);
%     ra = RMA(dataT,within,{0.05,{'bonferroni','hsd'}});
    ra = RMA(dataT,within,{0.05,''});
    ra.ranova
    print_for_manuscript(ra)
    %%
    redF = [2]; redV = {1};
    [dataTR,withinR] = reduce_within_between(dataT,within,redF,redV);
    raR = RMA(dataTR,withinR,{0.025,{'hsd'}});
    raR.ranova
    print_for_manuscript(raR)
    %%
    an = 4;
    cn = 3; tn = 4;
    respPVB_t = allresp_trialsC{an,1}; respPVB_i = allresp_trialsC{an,2};respPVB1 = cell_list_op(cell_list_op(respPVB_t,respPVB_i,'or'),[],'or',1);
    respPVB_t1 = allresp_trialsC{an,9}; respPVB_i1 = allresp_trialsC{an,10};respPVB2 = cell_list_op(cell_list_op(respPVB_t1,respPVB_i1,'or'),[],'or',1);
    
    respPVNB_t = allresp_trialsC{an,3}; respPVNB_i = allresp_trialsC{an,4}; respPVNB1 = cell_list_op(cell_list_op(respPVNB_t,respPVNB_i,'or'),[],'or',1);
    respPVNB_t1 = allresp_trialsC{an,5}; respPVNB_i1 = allresp_trialsC{an,6}; respPVNB2 = cell_list_op(cell_list_op(respPVNB_t1,respPVNB_i1,'or'),[],'or',1);
    respPVNB_t2 = allresp_trialsC{an,7}; respPVNB_i2 = allresp_trialsC{an,8}; respPVNB3 = cell_list_op(cell_list_op(respPVNB_t2,respPVNB_i2,'or'),[],'or',1);
    
    respPVB = cell_list_op(respPVB1,respPVB2,'or');
    respPVNB = cell_list_op(respPVNB1,respPVNB2,'or'); respPVNB = cell_list_op(respPVNB,respPVNB3,'or');
    
%     respPVB = cell_list_op(respPVB_t,respPVB_t1,'or'); respPVB = cell_list_op(respPVB,[],'or',1)
%     respPVNB = cell_list_op(respPVNB_t,respPVNB_t1,'or'); respPVNB = cell_list_op(respPVNB,respPVNB_t2,'or');    respPVNB = cell_list_op(respPVNB,[],'or',1)
    
    respPVBo = cell_list_op(respPVB,cell_list_op(respPVNB,[],'not'),'and');
    respPVNBo = cell_list_op(cell_list_op(respPVB,[],'not'),respPVNB,'and');
    respPVBNB = cell_list_op(respPVB,respPVNB,'and');
    
    [sum(respPVBo{1}) sum(respPVNBo{1}) sum(respPVBNB{1})]
    
    mR = allmRsT{cn}(an,:); respPV = allresp_trialsC{an,cn};
    respPVor = cell_list_op(respPV,[],'or',1);
%     respT = respPVor;
    respT = respPVBNB;
    [pv,pvc,pvct,pvg] = get_population_vector(mR,respT);
    [pva,cnr,pkpos] = align_by_peaks(pv{tn});
    [pva_all,cns,pkpos] = align_by_peaks(pv,cnr);
    
   
    
    ff = makeFigureRowsCols(108,[1.5 7 14 2.5],'RowsCols',[2 10],...
        'spaceRowsCols',[0.05 0.025],'rightUpShifts',[0.03 0.13],'widthHeightAdjustment',...
        [-30 -150]);  
    for tn = 1:10
%     hf = get_figure(8,[5 1 9 9]);hold on;
        pva = pva_all{tn};
%         [pva,cnr,pkpos] = align_by_peaks(pv{tn});
        pvab = zeros(size(pva)); pvab(pva==1) = 1; [yy,xx] = find(pvab);
        axes(ff.h_axes(1,tn));
        imagesc(pva);colorbar; set(gca,'Ydir','Normal');
        R = corrcoef(xx,yy);
        title(R(1,2));
        axes(ff.h_axes(2,tn));
        scatter(xx,yy);
        
    end