function processing_trial_by_trial
%%
% G = 'C'; o = oC;
% G = 'A'; o = oA;

si = [C1_t_D C1_i_T C2_t_D C2_i_T C3_t_D C3_i_T C4_t_D C4_i_T];
si = [C1_t_D C2_t_D C3_t_D C4_t_D];
[allRsC_C,allmRsT_C,allresp_C,pcs_C,event_type] = return_values_trial_to_trial_Analysis(oC,si);

[allRsC_A,allmRsT_A,allresp_A,pcs_A,~] = return_values_trial_to_trial_Analysis(oA,si);
%%
type = 3;
ff = makeFigureRowsCols(107,[1 0.5 4 1],'RowsCols',[1 5],'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.07 0.26],'widthHeightAdjustment',[10 -410]);
set(gcf,'color','w');    set(gcf,'Position',[10 3 6.9 1.9]);
MY = 8; ysp = 1; mY = 0; % responsive cells
stp = 0.11; widths = [1.3 1.3 1.3 1.3 1.3 0.5 0.5 0.5]-0.1; gap = 0.16;
adjust_axes(ff,[mY MY],stp,widths,gap,{'Cell #'});
% heatmap_conj_comp(hf,allresp_C(5,:),1,{si,rasterNamesTxt,'C'});
for ii = 1:length(ff.h_axes)
    heatmap_conj_comp(ff.h_axes(1,ii),allresp_C(ii,:),type,{si,rasterNamesTxt});
end
save_pdf(ff.hf,mData.pdf_folder,sprintf('OI_Map_mean_tu_spatial_C.pdf'),600);

ff = makeFigureRowsCols(107,[1 0.5 4 1],'RowsCols',[1 5],'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.07 0.26],'widthHeightAdjustment',[10 -410]);
set(gcf,'color','w');    set(gcf,'Position',[10 3 6.9 1.9]);
MY = 8; ysp = 1; mY = 0; % responsive cells
stp = 0.11; widths = [1.3 1.3 1.3 1.3 1.3 0.5 0.5 0.5]-0.1; gap = 0.16;
adjust_axes(ff,[mY MY],stp,widths,gap,{'Cell #'});
% heatmap_conj_comp(hf,allresp_C(5,:),1,{si,rasterNamesTxt,'C'});
for ii = 1:length(ff.h_axes)
    heatmap_conj_comp(ff.h_axes(1,ii),allresp_A(ii,:),type,{si,rasterNamesTxt});
end
save_pdf(ff.hf,mData.pdf_folder,sprintf('OI_Map_mean_tu_spatial_A.pdf'),600);


%% for responsivity
varC = pcs_C.respRV;
varA = pcs_A.respRV;
[within,dvn,xlabels] = make_within_table({'Cond','DT','Trials'},[4,2,10]);
dataT = make_between_table({varC;varA},dvn);
ra = RMA(dataT,within,{0.05,{'hsd','bonferroni'}});
ra.ranova


%% for conjunction

varC = []; varA = [];
for an = 1:5
    tvar = pcs_C.conjV(an,:);
%     tvar = pcs_C.comp1V(an,:);
%     tvar = pcs_C.comp2V(an,:);
    tvar(10:10:79) = NaN; tvar(isnan(tvar)) = []; 
    varC(an,:) = tvar;
    tvar = pcs_A.conjV(an,:);
%     tvar = pcs_A.comp1V(an,:);
%     tvar = pcs_A.comp2V(an,:);
    tvar(10:10:79) = NaN; tvar(isnan(tvar)) = []; 
    varA(an,:) = tvar;
end

[within,dvn,xlabels,withinD] = make_within_table({'Cond','DT','Trials'},[4,2,9]); awithinD = withinD;

    trialN = 1;
    varC2 = varC(:,awithinD(:,1) == 2 & awithinD(:,2) == 1);
    varA2 = varA(:,awithinD(:,1) == 2 & awithinD(:,2) == 1);
    
    [within1,dvn,xlabels] = make_within_table({'TrialPs'},[9]);
    dataT = make_between_table({varC2;varA2},dvn);
    ra = RMA(dataT,within1,{0.05,{'hsd','bonferroni'}});
    ra.ranova
%%
dataT = make_between_table({varC;varA},dvn);
ra = RMA(dataT,within,{0.05,{'hsd','bonferroni'}});
ra.ranova

%%
prop_names = {'resp','N_Resp_Trials','zMI','zMINaN','HaFD','HiFD','cells_pooled'};
event_type = {'1-D','2-D','3-D','4-D','1-T','2-T','3-T','4-T'};
sic = {[C1_t_D];[C2_t_D];[C3_t_D];[C4_t_D];[C1_i_T];[C2_i_T];[C3_i_T];[C4_i_T]};
pni = 7;
[all_gFR_C,all_gV_C] = return_values_props(oC,sic,pni);
[all_gFR_A,all_gV_A] = return_values_props(oA,sic,pni);

%%
while 1
    [~,varC] = find_trial_responsivity(allmRsT_C);
    varC1 = find_percent(allresp_C);
    varA1 = find_percent(allresp_A);
    
    [within,dvn,xlabels,withinD] = make_within_table({'Cond','DT','Trials'},[4,2,10]); awithinD = withinD;
    %%
    varC2 = varC1(:,awithinD(:,1) == 1 & awithinD(:,2) == 2);
    varA2 = varA1(:,awithinD(:,1) == 1 & awithinD(:,2) == 2);
    
    [within1,dvn,xlabels] = make_within_table({'Trials'},[10]);
    dataT = make_between_table({varC2;varA2},dvn);
    ra = RMA(dataT,within1,{0.05,{'hsd','bonferroni'}});
    ra.ranova
    
    %%
    varC2 = varC1(:,awithinD(:,2) == 2);
    varA2 = varA1(:,awithinD(:,2) == 2);
    
    [within1,dvn,xlabels] = make_within_table({'Cond','Trials'},[4 10]);
    dataT = make_between_table({varC2;varA2},dvn);
    ra = RMA(dataT,within1,{0.05,{'hsd','bonferroni'}});
    ra.ranova
    
     %%
     trialN = 1;
    varC2 = varC1(:,awithinD(:,2) == 1 & awithinD(:,3) == trialN);
    varA2 = varA1(:,awithinD(:,2) == 1 & awithinD(:,3) == trialN);
    
    [within1,dvn,xlabels] = make_within_table({'Cond'},[4]);
    dataT = make_between_table({varC2;varA2},dvn);
    ra = RMA(dataT,within1,{0.05,{'hsd','bonferroni'}});
    ra.ranova
    
    
    %%
    break;
end

%%
types = {'Conj','Comp1','Comp2'};
for ti = 1:length(types)
    hf = get_figure(6,[8 3 3.5 3.5]);
    mOI = heatmap_conj_comp(hf,allresp_A(:,:),ti,{si,rasterNamesTxt});
    save_pdf(hf,mData.pdf_folder,sprintf('OI_Map_mean_tu_spatial_A_%s.pdf',types{ti}),600);
    cluster_conj_comp(hf,mOI,event_type,[]);
    save_pdf(hf,mData.pdf_folder,sprintf('heatmap_cluster_A_%s.pdf',types{ti}),600);
end

for ti = 1:length(types)
    hf = get_figure(6,[8 3 3.5 3.5]);
    mOI = heatmap_conj_comp(hf,allresp_C(:,:),ti,{si,rasterNamesTxt});
    save_pdf(hf,mData.pdf_folder,sprintf('OI_Map_mean_tu_spatial_C_%s.pdf',types{ti}),600);
    cluster_conj_comp(hf,mOI,event_type,[]);
    save_pdf(hf,mData.pdf_folder,sprintf('heatmap_cluster_C_%s.pdf',types{ti}),600);
end