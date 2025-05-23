
%% population vector and correlation Distance 
while 1
    titles = {'Dur','Dis','Dur','Dis'};
    an = 4; cn = 1;
    cell_vars = {dur_cells_T dis_cells_T dur_cells_I dis_cells_I};
    si = [Ars_t_T Ars_t_D Ars_i_T Ars_i_D];
    si = [Ars_t_D Ars_t_D Ars_i_T Ars_i_T];
    Rs = o.Rs(:,si); mR = o.mR(:,si);
    props1 = get_props_Rs(Rs,50);
    cii = 4; cellvar = cell_vars{cii};
    resp = repmat(cellvar(:,cn),1,4);
    p_resp = find_percent(resp);

    ff = makeFigureRowsCols(107,[1 0.5 4 1],'RowsCols',[2 4],...
        'spaceRowsCols',[0 0.03],'rightUpShifts',[0.11 0.11],'widthHeightAdjustment',...
        [-55 -80]);    set(gcf,'color','w');    set(gcf,'Position',[10 3 3.5 1.5]);
    [CRc,aCRc,mRR] = find_population_vector_corr(Rs,mR,resp,cii);
    ff = show_population_vector_and_corr(mData,ff,Rs(an,:),mRR(an,:),CRc(an,:),[],[]);
    for ii = 1:length(ff.h_axes(1,:)) ht = get_obj(ff.h_axes(1,ii),'title'); set_obj(ht,{'String',titles{ii}}); end
    changePosition(ff.h_axes(2,1).YLabel,[-1 0 0]); 
    for ii = 1:length(ff.h_axes(1,:)) set_obj(ff.h_axes(2,ii),{'xtick',[1 68 136],'xticklabels',[0 7.5 15]}); set_obj(ff.h_axes(2,ii),{'ytick',[1 68 136],'yticklabels',[0 7.5 15]}); end 
    save_pdf(ff.hf,mData.pdf_folder,sprintf('PV_temporal.pdf'),600);

   % average correlation of all animals
    ff = makeFigureRowsCols(108,[1 0.5 4 0.5],'RowsCols',[1 4],...
        'spaceRowsCols',[0 0.03],'rightUpShifts',[0.11 0.25],'widthHeightAdjustment',...
        [-55 -350]);    set(gcf,'color','w');    set(gcf,'Position',[8 8 3.5 0.8]);
    ff = show_population_vector_and_corr(mData,ff,Rs(an,:),[],aCRc,[],[]);
    changePosition(ff.h_axes(1,1).YLabel,[-1 0 0]); 
    for ii = 1:length(ff.h_axes(1,:)) set_obj(ff.h_axes(1,ii),{'xtick',[1 68 136],'xticklabels',[0 7.5 15]}); set_obj(ff.h_axes(1,ii),{'ytick',[1 68 136],'yticklabels',[0 7.5 15]}); end 
    save_pdf(ff.hf,mData.pdf_folder,sprintf('aPV_temporal.pdf'),600);
    %%
    break;
end

%% population vector and correlation Distance 
while 1
    titles = {'Dur','Dis','Dur','Dis'};
    an = 4; cn = 2;
    cell_vars = {dur_cells_T dis_cells_T dur_cells_I dis_cells_I};
    si = [Ars_t_T Ars_t_D Ars_i_T Ars_i_D];
%     si = [Ars_t_D Ars_t_D Ars_i_T Ars_i_T];
    Rs = o.Rs(:,si); mR = o.mR(:,si);
    props1 = get_props_Rs(Rs,50);
    resp = [];
    for cii = 1:4
        cellvar = cell_vars{cii};
        resp = [resp cellvar(:,cn)];
    end
    p_resp = find_percent(resp);

    ff = makeFigureRowsCols(107,[1 0.5 4 1],'RowsCols',[2 4],...
        'spaceRowsCols',[0 0.03],'rightUpShifts',[0.11 0.11],'widthHeightAdjustment',...
        [-55 -80]);    set(gcf,'color','w');    set(gcf,'Position',[10 3 3.5 1.5]);
    [CRc,aCRc,mRR] = find_population_vector_corr(Rs,mR,resp,0);
    ff = show_population_vector_and_corr(mData,ff,Rs(an,:),mRR(an,:),CRc(an,:),[],[]);
    for ii = 1:length(ff.h_axes(1,:)) ht = get_obj(ff.h_axes(1,ii),'title'); set_obj(ht,{'String',titles{ii}}); end
    changePosition(ff.h_axes(2,1).YLabel,[-1 0 0]); 
    for ii = 1:length(ff.h_axes(1,:)) set_obj(ff.h_axes(2,ii),{'xtick',[1 68 136],'xticklabels',[0 7.5 15]}); set_obj(ff.h_axes(2,ii),{'ytick',[1 68 136],'yticklabels',[0 7.5 15]}); end 
    save_pdf(ff.hf,mData.pdf_folder,sprintf('PV_temporal.pdf'),600);

   % average correlation of all animals
    ff = makeFigureRowsCols(108,[1 0.5 4 0.5],'RowsCols',[1 4],...
        'spaceRowsCols',[0 0.03],'rightUpShifts',[0.11 0.25],'widthHeightAdjustment',...
        [-55 -350]);    set(gcf,'color','w');    set(gcf,'Position',[8 8 3.5 0.8]);
    ff = show_population_vector_and_corr(mData,ff,Rs(an,:),[],aCRc,[],[]);
    changePosition(ff.h_axes(1,1).YLabel,[-1 0 0]); 
    for ii = 1:length(ff.h_axes(1,:)) set_obj(ff.h_axes(1,ii),{'xtick',[1 68 136],'xticklabels',[0 7.5 15]}); set_obj(ff.h_axes(1,ii),{'ytick',[1 68 136],'yticklabels',[0 7.5 15]}); end 
    save_pdf(ff.hf,mData.pdf_folder,sprintf('aPV_temporal.pdf'),600);
    %%
    break;
end


%% population vector and correlation temporal 
while 1
    selected_property = 'tuned';
    titles = {'3','4','5'};
    an = 4;
%     si = [Ar_i_T ArL_i_T Ars_i_T];
    si = [Ar_On ArL_On Ars_On];
%     si = [Ar_Off ArL_Off Ars_Off];
    Rs = o.Rs(:,si);mR = o.mR(:,si);
    ntrials = 50;
    props1 = get_props_Rs(Rs,ntrials);
    resp = props1.vals;
    ff = makeFigureRowsCols(107,[1 0.5 4 1],'RowsCols',[2 3],...
        'spaceRowsCols',[0 0.03],'rightUpShifts',[0.11 0.11],'widthHeightAdjustment',...
        [-55 -80]);    set(gcf,'color','w');    set(gcf,'Position',[10 3 2.5 1.5]);
    [CRc,aCRc,mRR] = find_population_vector_corr(Rs,mR,resp,0);
    ff = show_population_vector_and_corr(mData,ff,Rs(an,:),mRR(an,:),CRc(an,:),[],[]);
%     for ii = 1:length(ff.h_axes(1,:)) ht = get_obj(ff.h_axes(1,ii),'title'); set_obj(ht,{'String',titles{ii}}); end
%     changePosition(ff.h_axes(2,1).YLabel,[-3 0 0]); 
%     for ii = 1:length(ff.h_axes(1,:)) set_obj(ff.h_axes(2,ii),{'xtick',[1 68 136],'xticklabels',[0 7.5 15]}); set_obj(ff.h_axes(2,ii),{'ytick',[1 68 136],'yticklabels',[0 7.5 15]}); end 
    save_pdf(ff.hf,mData.pdf_folder,sprintf('PV_temporal_%s.pdf',selected_property),600);

   % average correlation of all animals
    ff = makeFigureRowsCols(108,[1 0.5 4 0.5],'RowsCols',[1 3],...
        'spaceRowsCols',[0 0.03],'rightUpShifts',[0.11 0.25],'widthHeightAdjustment',...
        [-55 -350]);    set(gcf,'color','w');    set(gcf,'Position',[8 6 2.5 0.8]);
    ff = show_population_vector_and_corr(mData,ff,Rs(an,:),[],aCRc,[],[]);
    changePosition(ff.h_axes(1,1).YLabel,[-3 0 0]); 
    for ii = 1:length(ff.h_axes(1,:)) set_obj(ff.h_axes(1,ii),{'xtick',[1 68 136],'xticklabels',[0 7.5 15]}); set_obj(ff.h_axes(1,ii),{'ytick',[1 68 136],'yticklabels',[0 7.5 15]}); end 
    save_pdf(ff.hf,mData.pdf_folder,sprintf('aPV_temporal_%s.pdf',selected_property),600);
    %%
    break;
end


%% population vector and correlation spatial for presentation
while 1
    selected_property = 'tuned';
    titles = {'Ar-t-D','Ar-t-D','Ar-t-D'};
    an = 4;
    si = [Ar_t_D Ar_t_D Ar_t_D];
    Rs = o.Rs(:,si);mR = o.mR(:,si);
    ntrials = 50;
    props1 = get_props_Rs(Rs,ntrials);
    q_type = 'All';
    resp = [props1.good_FR(:,1) props1.good_zMI(:,1) props1.good_Gauss(:,1)]; 
%     resp = cell_list_op(props1.good_FR_and_Gauss_loose,[],'or');
%     resp = cell_list_op(props1.good_FR,[],'or');
    ff = makeFigureRowsCols(107,[1 0.5 4 1],'RowsCols',[2 3],...
        'spaceRowsCols',[0 0.03],'rightUpShifts',[0.11 0.11],'widthHeightAdjustment',...
        [-55 -80]);    set(gcf,'color','w');    set(gcf,'Position',[10 3 2.5 1.5]);
    [CRc,aCRc,mRR] = find_population_vector_corr(Rs,mR,resp,0);
    ff = show_population_vector_and_corr(mData,ff,Rs(an,:),mRR(an,:),CRc(an,:),[],[]);
    for ii = 1:length(ff.h_axes(1,:)) ht = get_obj(ff.h_axes(1,ii),'title'); set_obj(ht,{'String',titles{ii}}); end
    changePosition(ff.h_axes(2,1).YLabel,[-3 0 0]); 
    colormap parula
    save_pdf(ff.hf,mData.pdf_folder,sprintf('PV_spatial_%s_%s.pdf',selected_property,q_type),600);

   % average correlation of all animals
    ff = makeFigureRowsCols(108,[1 0.5 4 0.5],'RowsCols',[1 3],...
        'spaceRowsCols',[0 0.03],'rightUpShifts',[0.11 0.25],'widthHeightAdjustment',...
        [-55 -350]);    set(gcf,'color','w');    set(gcf,'Position',[8 8 2.5 0.8]);
    ff = show_population_vector_and_corr(mData,ff,Rs(an,:),[],aCRc,[],[]);
    changePosition(ff.h_axes(1,1).YLabel,[-3 0 0]); 
    colormap parula
    save_pdf(ff.hf,mData.pdf_folder,sprintf('aPV_spatial_%s_%s.pdf',selected_property,q_type),600);
    %%
    break;
end


%% population vector and correlation spatial 
while 1
    selected_property = 'tuned';
    titles = {'3 ','4 ','5 '};
    an = 4;
    si = [Ar_t_D ArL_t_D Ars_t_D];
    si = [Ar_i_T ArL_i_T Ars_i_T];
    Rs = o.Rs(:,si);mR = o.mR(:,si);
    props1 = get_props_Rs(Rs,[50,100]);
    q_type = '1040';
    resp = props1.vals;
    ff = makeFigureRowsCols(107,[1 0.5 4 1],'RowsCols',[2 3],...
        'spaceRowsCols',[0 0.03],'rightUpShifts',[0.11 0.11],'widthHeightAdjustment',...
        [-55 -80]);    set(gcf,'color','w');    set(gcf,'Position',[10 3 2.5 1.5]);
    [CRc,aCRc,mRR] = find_population_vector_corr(Rs,mR,resp,0);
    ff = show_population_vector_and_corr(mData,ff,Rs(an,:),mRR(an,:),CRc(an,:),[],[]);
    for ii = 1:length(ff.h_axes(1,:)) ht = get_obj(ff.h_axes(1,ii),'title'); set_obj(ht,{'String',titles{ii}}); end
    axes(ff.h_axes(1,2));hold on;
    ylims = ylim;
%     plot([36,36],ylims,'m');
    
    changePosition(ff.h_axes(2,1).YLabel,[-3 0 0]);
    colormap_ig
    save_pdf(ff.hf,mData.pdf_folder,sprintf('PV_spatial_%s_%s.pdf',selected_property,q_type),600);

   % average correlation of all animals
    ff = makeFigureRowsCols(108,[1 0.5 4 0.5],'RowsCols',[1 3],...
        'spaceRowsCols',[0 0.03],'rightUpShifts',[0.11 0.25],'widthHeightAdjustment',...
        [-55 -350]);    set(gcf,'color','w');    set(gcf,'Position',[8 8 2.5 0.85]);
    ff = show_population_vector_and_corr(mData,ff,Rs(an,:),[],aCRc,[],[]);
    changePosition(ff.h_axes(1,1).YLabel,[-3 0 0]); 
    colormap_ig
    save_pdf(ff.hf,mData.pdf_folder,sprintf('aPV_spatial_%s_%s.pdf',selected_property,q_type),600);
    %%
    break;
end
