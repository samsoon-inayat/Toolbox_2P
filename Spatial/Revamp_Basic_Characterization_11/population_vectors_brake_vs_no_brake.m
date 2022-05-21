%% population vector and correlation sensory air
while 1
    selected_property = 'untuned';
    an = 1;
    titles = {'Air','Light','NB-A','NB-AL','NB-A*'};
    si = [Ab_On Lb Abs_On Ar_On ArL_On Ars_On];
%     si = [Ab_Off Abs_Off Ar_Off ArL_Off Ars_Off];
%     si = [Ar_D ArL_D Ars_D Ar_T ArL_T];
    Rs = o.Rs(:,si);mR = o.mR(:,si);
    ntrials = 50;
    props1 = get_props_Rs(Rs,ntrials);
    respA = props1.good_FR_and_untuned;
%     respA = props1.good_zMI;
    resp1 = cell_list_op(respA(:,1:2),[],'or');
    resp2 = cell_list_op(respA(:,3:4),[],'or');
    resp = respA;%[resp1(:,1:2),resp2(:,1:2)];
%     eval(cmdTxt);
    ffM = makeFigureRowsCols(107,[1 0.5 4 1],'RowsCols',[2 5],...
        'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.08 0.11],'widthHeightAdjustment',...
        [10 -80]);    set(gcf,'color','w');    set(gcf,'Position',[10 3 5 2]);
    sInds = 1:5;
    [CRc,aCRc1,mRR] = find_population_vector_corr(Rs(:,sInds),mR(:,sInds),resp(:,sInds),0);
    ff = ffM; ff.axesPos = ffM.axesPos(:,sInds); ff.h_axes = ffM.h_axes(:,sInds);
    ff = show_population_vector_and_corr(mData,ff,Rs(an,sInds),mRR(an,:),CRc(an,:),[],[],0);
    for ii = 1:length(ff.h_axes(1,:)) ht = get_obj(ff.h_axes(1,ii),'title'); set(ht,'String',titles{ii}); end
 
    save_pdf(ff.hf,mData.pdf_folder,sprintf('PV_air_%s.pdf',selected_property),600);

    % average correlation of all animals
    ffM = makeFigureRowsCols(108,[1 0.5 4 1],'RowsCols',[1 5],...
        'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.08 0.24],'widthHeightAdjustment',...
        [10 -300]);    set(gcf,'color','w');    set(gcf,'Position',[10 7 5 1]);
    sInds = 1:5;
    ff = ffM; ff.axesPos = ffM.axesPos(:,sInds); ff.h_axes = ffM.h_axes(:,sInds);
    ff = show_population_vector_and_corr(mData,ff,Rs(an,sInds),[],aCRc1,[],[],0);
%     sInds = 3:4;
%     ff = ffM; ff.axesPos = ffM.axesPos(:,sInds); ff.h_axes = ffM.h_axes(:,sInds);
%     ff = show_population_vector_and_corr(mData,ff,Rs(an,sInds),[],aCRc1,[],[]);
    save_pdf(ff.hf,mData.pdf_folder,sprintf('aPV_air_%s.pdf',selected_property),600);
    %%
    break;
end

%% population vector and correlation sensory light
while 1
    selected_property = 'untuned';
    an = 1;
    titles = {'B-L','B-L*','NB-L'};
    si = [Lb Lbs ArL_L];
    Rs = o.Rs(:,si);mR = o.mR(:,si);
    ntrials = 50;
    props1 = get_props_Rs(Rs,ntrials);
    respA = props1.good_FR_and_untuned;
    resp = respA;%[resp1(:,1:2),resp2(:,1:2)];
%     eval(cmdTxt);
    ffM = makeFigureRowsCols(107,[1 0.5 4 1],'RowsCols',[2 3],...
        'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.08 0.11],'widthHeightAdjustment',...
        [10 -80]);    set(gcf,'color','w');    set(gcf,'Position',[10 3 2.75 2]);
    sInds = 1:3;
    [CRc,aCRc1,mRR] = find_population_vector_corr(Rs(:,sInds),mR(:,sInds),resp(:,sInds),0);
    ff = ffM; ff.axesPos = ffM.axesPos(:,sInds); ff.h_axes = ffM.h_axes(:,sInds);
    ff = show_population_vector_and_corr(mData,ff,Rs(an,sInds),mRR(an,:),CRc(an,:),[],[],0);
    for ii = 1:length(ff.h_axes(1,:)) ht = get_obj(ff.h_axes(1,ii),'title'); set(ht,'String',titles{ii}); end
    
%     sInds = 3:4;
%     [CRc,aCRc2,mRR] = find_population_vector_corr(Rs(:,sInds),mR(:,sInds),resp(:,sInds),0);
%     ff = ffM; ff.axesPos = ffM.axesPos(:,sInds); ff.h_axes = ffM.h_axes(:,sInds);
%     ff = show_population_vector_and_corr(mData,ff,Rs(an,sInds),mRR(an,:),CRc(an,:),[],[]);
%     for ii = 1:length(ff.h_axes(1,:)) ht = get_obj(ff.h_axes(1,ii),'title'); set(ht,'String',titles{ii}); end
    
    
    save_pdf(ff.hf,mData.pdf_folder,sprintf('PV_air_%s.pdf',selected_property),600);

    % average correlation of all animals
    ffM = makeFigureRowsCols(108,[1 0.5 4 1],'RowsCols',[1 3],...
        'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.08 0.24],'widthHeightAdjustment',...
        [10 -300]);    set(gcf,'color','w');    set(gcf,'Position',[10 7 2.75 1]);
    sInds = 1:3;
    ff = ffM; ff.axesPos = ffM.axesPos(:,sInds); ff.h_axes = ffM.h_axes(:,sInds);
    ff = show_population_vector_and_corr(mData,ff,Rs(an,sInds),[],aCRc1,[],[],0);

    save_pdf(ff.hf,mData.pdf_folder,sprintf('aPV_air_%s.pdf',selected_property),600);
    %%
    break;
end


%% population vector and correlation sensory light + Control light
while 1
    selected_property = 'untuned';
    an = 5;
    titles = {'B-L','B-L*','NB-L'};
    si = [Ar_L ArL_L Ars_L];
    Rs = o.Rs(:,si);mR = o.mR(:,si);
    ntrials = 50;
    props1 = get_props_Rs(Rs,ntrials);
    respA = props1.vals;
    resp = respA;%[resp1(:,1:2),resp2(:,1:2)];
%     eval(cmdTxt);
    ffM = makeFigureRowsCols(107,[1 0.5 4 1],'RowsCols',[2 3],...
        'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.08 0.11],'widthHeightAdjustment',...
        [10 -80]);    set(gcf,'color','w');    set(gcf,'Position',[10 3 2.75 2]);
    sInds = 1:3;
    [CRc,aCRc1,mRR] = find_population_vector_corr(Rs(:,sInds),mR(:,sInds),resp(:,sInds),0);
    ff = ffM; ff.axesPos = ffM.axesPos(:,sInds); ff.h_axes = ffM.h_axes(:,sInds);
    ff = show_population_vector_and_corr(mData,ff,Rs(an,sInds),mRR(an,:),CRc(an,:),[],[],0);
    for ii = 1:length(ff.h_axes(1,:)) ht = get_obj(ff.h_axes(1,ii),'title'); set(ht,'String',titles{ii}); end
    
%     sInds = 3:4;
%     [CRc,aCRc2,mRR] = find_population_vector_corr(Rs(:,sInds),mR(:,sInds),resp(:,sInds),0);
%     ff = ffM; ff.axesPos = ffM.axesPos(:,sInds); ff.h_axes = ffM.h_axes(:,sInds);
%     ff = show_population_vector_and_corr(mData,ff,Rs(an,sInds),mRR(an,:),CRc(an,:),[],[]);
%     for ii = 1:length(ff.h_axes(1,:)) ht = get_obj(ff.h_axes(1,ii),'title'); set(ht,'String',titles{ii}); end
    
    
    save_pdf(ff.hf,mData.pdf_folder,sprintf('PV_air_%s.pdf',selected_property),600);

    % average correlation of all animals
    ffM = makeFigureRowsCols(108,[1 0.5 4 1],'RowsCols',[1 3],...
        'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.08 0.24],'widthHeightAdjustment',...
        [10 -300]);    set(gcf,'color','w');    set(gcf,'Position',[10 7 2.75 1]);
    sInds = 1:3;
    ff = ffM; ff.axesPos = ffM.axesPos(:,sInds); ff.h_axes = ffM.h_axes(:,sInds);
    ff = show_population_vector_and_corr(mData,ff,Rs(an,sInds),[],aCRc1,[],[],0);

    save_pdf(ff.hf,mData.pdf_folder,sprintf('aPV_air_%s.pdf',selected_property),600);
    %%
    break;
end



%% population vector and correlation distance
while 1
    selected_property = 'untuned';
    an = 1;
     titles = {'NB-A','NB-AL','NB-A*'};
    si = [Ar_D ArL_D Ars_D];
    Rs = o.Rs(:,si);mR = o.mR(:,si);
    ntrials = 50;
    props1 = get_props_Rs(Rs,ntrials);
    resp = props1.good_zMI;%[resp1(:,1:2),resp2(:,1:2)];
%     eval(cmdTxt);
    ffM = makeFigureRowsCols(107,[1 0.5 4 1],'RowsCols',[2 3],...
        'spaceRowsCols',[0.01 0.01],'rightUpShifts',[0.08 0.11],'widthHeightAdjustment',...
        [-20 -80]);    set(gcf,'color','w');    set(gcf,'Position',[10 3 2.75 2]);
    sInds = 1:3;
    [CRc,aCRc1,mRR] = find_population_vector_corr(Rs(:,sInds),mR(:,sInds),resp(:,sInds),0);
    ff = ffM; ff.axesPos = ffM.axesPos(:,sInds); ff.h_axes = ffM.h_axes(:,sInds);
    ff = show_population_vector_and_corr(mData,ff,Rs(an,sInds),mRR(an,:),CRc(an,:),[],[],0);
    for ii = 1:length(ff.h_axes(1,:)) ht = get_obj(ff.h_axes(1,ii),'title'); set(ht,'String',titles{ii}); end
    
%     sInds = 3:4;
%     [CRc,aCRc2,mRR] = find_population_vector_corr(Rs(:,sInds),mR(:,sInds),resp(:,sInds),0);
%     ff = ffM; ff.axesPos = ffM.axesPos(:,sInds); ff.h_axes = ffM.h_axes(:,sInds);
%     ff = show_population_vector_and_corr(mData,ff,Rs(an,sInds),mRR(an,:),CRc(an,:),[],[]);
%     for ii = 1:length(ff.h_axes(1,:)) ht = get_obj(ff.h_axes(1,ii),'title'); set(ht,'String',titles{ii}); end
    
    
    save_pdf(ff.hf,mData.pdf_folder,sprintf('PV_air_%s.pdf',selected_property),600);

    % average correlation of all animals
    ffM = makeFigureRowsCols(108,[1 0.5 4 1],'RowsCols',[1 3],...
        'spaceRowsCols',[0.01 0.01],'rightUpShifts',[0.08 0.24],'widthHeightAdjustment',...
        [-20 -300]);    set(gcf,'color','w');    set(gcf,'Position',[10 7 2.75 1]);
    sInds = 1:3;
    ff = ffM; ff.axesPos = ffM.axesPos(:,sInds); ff.h_axes = ffM.h_axes(:,sInds);
    ff = show_population_vector_and_corr(mData,ff,Rs(an,sInds),[],aCRc1,[],[],0);

    save_pdf(ff.hf,mData.pdf_folder,sprintf('aPV_air_%s.pdf',selected_property),600);
    %%
    break;
end


%% population vector and correlation distance
while 1
    selected_property = 'untuned';
    an = 1;
     titles = {'NB-A','NB-AL','NB-A*'};
    si = [Ar_T ArL_T Ars_T];
    Rs = o.Rs(:,si);mR = o.mR(:,si);
    ntrials = 50;
    props1 = get_props_Rs(Rs,ntrials);
    resp = props1.good_zMI;%[resp1(:,1:2),resp2(:,1:2)];
%     eval(cmdTxt);
    ffM = makeFigureRowsCols(107,[1 0.5 4 1],'RowsCols',[2 3],...
        'spaceRowsCols',[0.01 0.01],'rightUpShifts',[0.08 0.11],'widthHeightAdjustment',...
        [-20 -80]);    set(gcf,'color','w');    set(gcf,'Position',[10 3 2.75 2]);
    sInds = 1:3;
    [CRc,aCRc1,mRR] = find_population_vector_corr(Rs(:,sInds),mR(:,sInds),resp(:,sInds),0);
    ff = ffM; ff.axesPos = ffM.axesPos(:,sInds); ff.h_axes = ffM.h_axes(:,sInds);
    ff = show_population_vector_and_corr(mData,ff,Rs(an,sInds),mRR(an,:),CRc(an,:),[],[],0);
    for ii = 1:length(ff.h_axes(1,:)) ht = get_obj(ff.h_axes(1,ii),'title'); set(ht,'String',titles{ii}); end
    
%     sInds = 3:4;
%     [CRc,aCRc2,mRR] = find_population_vector_corr(Rs(:,sInds),mR(:,sInds),resp(:,sInds),0);
%     ff = ffM; ff.axesPos = ffM.axesPos(:,sInds); ff.h_axes = ffM.h_axes(:,sInds);
%     ff = show_population_vector_and_corr(mData,ff,Rs(an,sInds),mRR(an,:),CRc(an,:),[],[]);
%     for ii = 1:length(ff.h_axes(1,:)) ht = get_obj(ff.h_axes(1,ii),'title'); set(ht,'String',titles{ii}); end
    
    
    save_pdf(ff.hf,mData.pdf_folder,sprintf('PV_air_%s.pdf',selected_property),600);

    % average correlation of all animals
    ffM = makeFigureRowsCols(108,[1 0.5 4 1],'RowsCols',[1 3],...
        'spaceRowsCols',[0.01 0.01],'rightUpShifts',[0.08 0.24],'widthHeightAdjustment',...
        [-20 -300]);    set(gcf,'color','w');    set(gcf,'Position',[10 7 2.75 1]);
    sInds = 1:3;
    ff = ffM; ff.axesPos = ffM.axesPos(:,sInds); ff.h_axes = ffM.h_axes(:,sInds);
    ff = show_population_vector_and_corr(mData,ff,Rs(an,sInds),[],aCRc1,[],[],0);

    save_pdf(ff.hf,mData.pdf_folder,sprintf('aPV_air_%s.pdf',selected_property),600);
    %%
    break;
end





%% population vector and correlation distance and time
while 1
    an = 1;
    titles = {'NB-A','NB-AL','NB-A*','NB-A','NB-AL','NB-A*'};
    si = [Ar_D ArL_D Ars_D Ar_T ArL_T Ars_T];
    Rs = o.Rs(:,si);mR = o.mR(:,si);
    ntrials = 50;
    props1 = get_props_Rs(Rs,ntrials);
    resp = props1.good_zMI;
%     resp = respA;%[resp1(:,1:2),resp2(:,1:2)];
%     resp = [respDT.exc respDT.exc];
%     resp = [respDT.inh respDT.inh];
%     eval(cmdTxt);
    ffM = makeFigureRowsCols(107,[1 0.5 4 1],'RowsCols',[2 6],...
        'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.08 0.11],'widthHeightAdjustment',...
        [10 -80]);    set(gcf,'color','w');    set(gcf,'Position',[10 3 6 2]);
    sInds = 1:6;
    [CRc,aCRc1,mRR] = find_population_vector_corr(Rs(:,sInds),mR(:,sInds),resp(:,sInds),0);
    ff = ffM; ff.axesPos = ffM.axesPos(:,sInds); ff.h_axes = ffM.h_axes(:,sInds);
    ff = show_population_vector_and_corr(mData,ff,Rs(an,sInds),mRR(an,:),CRc(an,:),[],[],0);
    for ii = 1:length(ff.h_axes(1,:)) ht = get_obj(ff.h_axes(1,ii),'title'); set(ht,'String',titles{ii}); end
    for ii = 4:6 set_obj(ff.h_axes(2,ii),{'xtick',[1 68 136],'xticklabels',[0 7.5 15]}); set_obj(ff.h_axes(2,ii),{'ytick',[1 68 136],'yticklabels',[0 7.5 15]}); end 
    save_pdf(ff.hf,mData.pdf_folder,sprintf('aPV_spatial_temp.pdf'),600);

    % average correlation of all animals
    ffM = makeFigureRowsCols(108,[1 0.5 4 1],'RowsCols',[1 6],...
        'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.08 0.24],'widthHeightAdjustment',...
        [10 -300]);    set(gcf,'color','w');    set(gcf,'Position',[10 7 6 1]);
    sInds = 1:6;
    ff = ffM; ff.axesPos = ffM.axesPos(:,sInds); ff.h_axes = ffM.h_axes(:,sInds);
    ff = show_population_vector_and_corr(mData,ff,Rs(an,sInds),[],aCRc1,[],[],0);
    for ii = 4:6 set_obj(ff.h_axes(1,ii),{'xtick',[1 68 136],'xticklabels',[0 7.5 15]}); set_obj(ff.h_axes(1,ii),{'ytick',[1 68 136],'yticklabels',[0 7.5 15]}); end 
    save_pdf(ff.hf,mData.pdf_folder,sprintf('aPV_spatial_temp.pdf'),600);
    %%
    break;
end
