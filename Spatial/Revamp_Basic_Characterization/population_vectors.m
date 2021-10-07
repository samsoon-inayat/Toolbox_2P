%% run conjunctive_representation file first to load data
selected_property = 'good_FR_and_Gauss';
cmdTxt = sprintf('good_FR = props1.%s;',selected_property);

%% population vector and correlation sensory
while 1
    an = 4;
    titles = {'L','ArL-L','L*'};
    si = si_seq([1 11 9]);
    Rs = o.Rs(:,si);mR = o.mR(:,si);
    ntrials = 50;
    props1 = get_props_Rs(Rs,ntrials);
    eval(cmdTxt);
    ff = makeFigureRowsCols(107,[1 0.5 4 1],'RowsCols',[2 3],...
        'spaceRowsCols',[0 -0.02],'rightUpShifts',[0.08 0.1],'widthHeightAdjustment',...
        [-20 -60]);    set(gcf,'color','w');    set(gcf,'Position',[10 3 3 2]);
    [CRc,aCRc,mRR] = find_population_vector_corr(Rs,mR,good_FR,0);
    ff = show_population_vector_and_corr(mData,ff,Rs(an,:),mRR(an,:),CRc(an,:),[],[]);
    for ii = 1:length(ff.h_axes(1,:)) ht = get_obj(ff.h_axes(1,ii),'title'); set(ht,'String',titles{ii}); end
    save_pdf(ff.hf,mData.pdf_folder,sprintf('PV_light_%s.pdf',selected_property),600);

    % average correlation of all animals
    ff = makeFigureRowsCols(108,[1 0.5 4 0.5],'RowsCols',[1 3],...
        'spaceRowsCols',[0 -0.02],'rightUpShifts',[0.08 0.2],'widthHeightAdjustment',...
        [-20 -240]);    set(gcf,'color','w');    set(gcf,'Position',[10 8 3 1]);
    ff = show_population_vector_and_corr(mData,ff,Rs(an,:),[],aCRc,[],[]);
    save_pdf(ff.hf,mData.pdf_folder,sprintf('aPV_light_%s.pdf',selected_property),600);
    %%
    break;
end

%% population vector and correlation sensory
while 1
    an = 4;
    titles = {'A','A*'};
    si = si_seq([2 10]);
    Rs = o.Rs(:,si);mR = o.mR(:,si);
    ntrials = 50;
    props1 = get_props_Rs(Rs,ntrials);
    eval(cmdTxt);
    ff = makeFigureRowsCols(107,[1 0.5 4 1],'RowsCols',[2 2],...
        'spaceRowsCols',[0 -0.02],'rightUpShifts',[0.2 0.1],'widthHeightAdjustment',...
        [-150 -60]);    set(gcf,'color','w');    set(gcf,'Position',[10 3 3 2]);
    [CRc,aCRc,mRR] = find_population_vector_corr(Rs,mR,good_FR,0);
    ff = show_population_vector_and_corr(mData,ff,Rs(an,:),mRR(an,:),CRc(an,:),[],[]);
    for ii = 1:length(ff.h_axes(1,:)) ht = get_obj(ff.h_axes(1,ii),'title'); set(ht,'String',titles{ii}); end
    save_pdf(ff.hf,mData.pdf_folder,sprintf('PV_air_%s.pdf',selected_property),600);

    % average correlation of all animals
    ff = makeFigureRowsCols(108,[1 0.5 4 0.5],'RowsCols',[1 2],...
        'spaceRowsCols',[0 -0.02],'rightUpShifts',[0.2 0.2],'widthHeightAdjustment',...
        [-150 -240]);    set(gcf,'color','w');    set(gcf,'Position',[10 8 3 1]);
    ff = show_population_vector_and_corr(mData,ff,Rs(an,:),[],aCRc,[],[]);
    save_pdf(ff.hf,mData.pdf_folder,sprintf('aPV_air_%s.pdf',selected_property),600);
    %%
    break;
end

%% population vector and correlation  temporal
while 1
    titles = {'Ar-it-Time','ArL-it-Time','Ar*-it-Time'};
    an = 4;
    si = si_seq(setdiff(1:11,[1 11 9 2 10]));
    si = si([2 4 6]);
    Rs = o.Rs(:,si);mR = o.mR(:,si);
    ntrials = 50;
    props1 = get_props_Rs(Rs,ntrials);
    eval(cmdTxt);
    ff = makeFigureRowsCols(107,[1 0.5 4 1],'RowsCols',[2 3],...
        'spaceRowsCols',[0 0.01],'rightUpShifts',[0.11 0.1],'widthHeightAdjustment',...
        [-35 -60]);    set(gcf,'color','w');    set(gcf,'Position',[10 3 3.45 2]);
    [CRc,aCRc,mRR] = find_population_vector_corr(Rs,mR,good_FR,0);
    ff = show_population_vector_and_corr(mData,ff,Rs(an,:),mRR(an,:),CRc(an,:),[],[]);
    for ii = 1:length(ff.h_axes(1,:)) ht = get_obj(ff.h_axes(1,ii),'title'); set_obj(ht,{'String',titles{ii}}); end
    save_pdf(ff.hf,mData.pdf_folder,sprintf('PV_temporal_%s.pdf',selected_property),600);

   % average correlation of all animals
    ff = makeFigureRowsCols(108,[1 0.5 4 0.5],'RowsCols',[1 3],...
        'spaceRowsCols',[0 0.01],'rightUpShifts',[0.11 0.2],'widthHeightAdjustment',...
        [-35 -240]);    set(gcf,'color','w');    set(gcf,'Position',[5 8 3.45 1]);
    ff = show_population_vector_and_corr(mData,ff,Rs(an,:),[],aCRc,[],[]);
    save_pdf(ff.hf,mData.pdf_folder,sprintf('aPV_temporal_%s.pdf',selected_property),600);
    %%
    break;
end


%% population vector and correlation spatial 
while 1
    titles = {'Ar-t-Dist','ArL-t-Dist','Ar*-t-Dist'};
    an = 4;
    si = si_seq(setdiff(1:11,[1 11 9 2 10]));
    si = si([1 3 5]);
    Rs = o.Rs(:,si);mR = o.mR(:,si);
    ntrials = 50;
    props1 = get_props_Rs(Rs,ntrials);
    eval(cmdTxt);
    ff = makeFigureRowsCols(107,[1 0.5 4 1],'RowsCols',[2 3],...
        'spaceRowsCols',[0 0.01],'rightUpShifts',[0.11 0.1],'widthHeightAdjustment',...
        [-35 -60]);    set(gcf,'color','w');    set(gcf,'Position',[10 3 3.45 2]);
    [CRc,aCRc,mRR] = find_population_vector_corr(Rs,mR,good_FR,0);
    ff = show_population_vector_and_corr(mData,ff,Rs(an,:),mRR(an,:),CRc(an,:),[],[]);
    for ii = 1:length(ff.h_axes(1,:)) ht = get_obj(ff.h_axes(1,ii),'title'); set_obj(ht,{'String',titles{ii}}); end
    save_pdf(ff.hf,mData.pdf_folder,sprintf('PV_spatial_%s.pdf',selected_property),600);

   % average correlation of all animals
    ff = makeFigureRowsCols(108,[1 0.5 4 0.5],'RowsCols',[1 3],...
        'spaceRowsCols',[0 0.01],'rightUpShifts',[0.11 0.2],'widthHeightAdjustment',...
        [-35 -240]);    set(gcf,'color','w');    set(gcf,'Position',[5 8 3.45 1]);
    ff = show_population_vector_and_corr(mData,ff,Rs(an,:),[],aCRc,[],[]);
    save_pdf(ff.hf,mData.pdf_folder,sprintf('aPV_spatial_%s.pdf',selected_property),600);
    %%
    break;
end

%% population vector and correlation spatial temporal
while 1
    titles = {'C3t-Dist','C4t-Dist','C3''t-Dist','C3it-Dist','C4it-Dist','C3''it-Dist'};
    an = 2;
    si1 = si_seq(setdiff(1:11,[1 11 9 2 10]));
    Rs1 = o.Rs(:,si1); mR1 = o.mR(:,si1);     ntrials = 50;     props1 = get_props_Rs(Rs1,ntrials);

    si = [si_air_dist_trials si_air_dist_itrials];
    Rs = o.Rs(:,si);mR = o.mR(:,si);
    eval(cmdTxt);
    ff = makeFigureRowsCols(107,[1 0.5 4 1],'RowsCols',[2 6],...
        'spaceRowsCols',[0 -0.03],'rightUpShifts',[0.05 0.1],'widthHeightAdjustment',...
        [25 -60]);    set(gcf,'color','w');    set(gcf,'Position',[10 3 6.99 2]);
    [CRc,aCRc,mRR] = find_population_vector_corr(Rs,mR,good_FR,0);
    ff = show_population_vector_and_corr(mData,ff,Rs(an,:),mRR(an,:),CRc(an,:),[],[]);
    for ii = 1:length(ff.h_axes(1,:)) ht = get_obj(ff.h_axes(1,ii),'title'); set_obj(ht,{'String',titles{ii}}); end
%     save_pdf(ff.hf,mData.pdf_folder,sprintf('population_vector_corr_time.pdf'),600);

   % average correlation of all animals
    ff = makeFigureRowsCols(108,[1 0.5 4 0.5],'RowsCols',[1 6],...
        'spaceRowsCols',[0 -0.03],'rightUpShifts',[0.05 0.2],'widthHeightAdjustment',...
        [20 -240]);    set(gcf,'color','w');    set(gcf,'Position',[5 8 6.99 1]);
    ff = show_population_vector_and_corr(mData,ff,Rs(an,:),[],aCRc,[],[]);
%     save_pdf(ff.hf,mData.pdf_folder,sprintf('average_population_vector_corr_time.pdf'),600);
    %%
    break;
end
