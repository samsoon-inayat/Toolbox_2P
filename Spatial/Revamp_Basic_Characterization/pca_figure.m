function pca_figure

mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size;
ei_11_15 = evalin('base','ei'); 

selContexts = [2 7];
rasterNames = {'air55T','air55T'};
RsR = get_rasters_data(ei_11_15,selContexts,rasterNames);
Rs = get_rasters_data_pca(ei_11_15,selContexts,rasterNames);
% Rs = get_rasters_data(ei_2_3,selContexts,rasterNames);
Rs = find_responsive_rasters(Rs,1:10);
mR = calc_mean_rasters(Rs,1:10);
[resp_fractionC,resp_valsC,OIC,mean_OIC,resp_ORC,resp_OR_fractionC,resp_ANDC,resp_AND_fractionC,resp_exc_inh] = get_responsive_fraction(Rs);

[respE,respE_OR,respE_AND,respE_fraction] = get_cell_list_exc_inh(resp_exc_inh,1,1);
[CRcE,aCRcE,mRRE] = find_population_vector_corr(Rs,mR,respE,0);

n = 0;
%% population vector and correlation single animal
if 1
    an = 1;
    ff = makeFigureRowsCols(106,[1 0.5 4 1],'RowsCols',[2 2],...
        'spaceRowsCols',[0 -0.02],'rightUpShifts',[0.15 0.1],'widthHeightAdjustment',...
        [-80 -70]);
    set(gcf,'color','w');
    set(gcf,'Position',[5 5 2.2 2]);
    [resp_fractionC,resp_valsC,OIC,mean_OIC,resp_ORC,resp_OR_fractionC,resp_ANDC,resp_AND_fractionC,resp_exc_inh] = get_responsive_fraction(Rs);
    resp = get_cell_list(resp_valsC,[1;2]);
    [CRc,aCRc,mRR] = find_population_vector_corr(Rs,mR,0,0);
    % ff = show_population_vector_and_corr(mData,ff,Rs(an,:),mRR(an,:),CRc(an,:),[-0.1 1],[]);
    ff = show_population_vector_and_corr(mData,ff,Rs(an,:),mRR(an,:),CRc(an,:),[-0.1 1],[]);
    set_obj(ff,{'FontWeight','Normal','FontSize',6,'LineWidth',0.25});
    ht = get_obj(ff,'title'); hyl = get_obj(ff,'ylabel'); changePosition(hyl(1,1),[7 0 0]);
    set_obj(ht,'String',{'Pop. Activity','Pop. Activity';'Pop. Correlation','Pop.Correlation'});
    set_obj(ht,{'FontSize',5,'FontWeight','Normal'});
    save_pdf(ff.hf,mData.pdf_folder,sprintf('pca_population_vector_corr.pdf'),600);
end
%% average population correlation (from all animals)
if 1
    ff = makeFigureRowsCols(107,[1 0.5 4 0.5],'RowsCols',[1 2],...
        'spaceRowsCols',[0 -0.03],'rightUpShifts',[0.15 0.2],'widthHeightAdjustment',...
        [-70 -240]);
    set(gcf,'color','w');
    set(gcf,'Position',[5 3 2.2 1]);
    ff = show_population_vector_and_corr(mData,ff,Rs(an,:),[],aCRc,[-0.1 1],[]);
    set_obj(ff,{'FontWeight','Normal','FontSize',6,'LineWidth',0.25});
    ht = get_obj(ff,'title'); hyl = get_obj(ff,'ylabel'); changePosition(hyl(1,1),[7 0 0]);
    set_obj(ht,'String',{'Avg. Pop. Correlation','Avg. Pop.Correlation'});
    set_obj(ht,{'FontSize',5,'FontWeight','Normal'});
    save_pdf(ff.hf,mData.pdf_folder,sprintf('pca_air_average_population_vector_corr.pdf'),600);
end