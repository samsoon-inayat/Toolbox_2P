function population_vector_corr

mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size;
ei_C = evalin('base','ei10_C2'); 
ei_A = evalin('base','ei10_A2'); 


selContexts = [1 2 3 4];
rasterNames = {'airD','airD','airD','airD'};

RsC = get_rasters_data(ei_C,selContexts,rasterNames);
mRsC = calc_mean_rasters(RsC,1:10);
RsC = find_responsive_rasters(RsC,1:10);
[CRcC,aCRcC,mRsRC] = find_population_vector_corr(RsC,mRsC,1);
% view_population_vector(Rs,mRs,300);
[resp_fractionC,resp_valsC,OIC,mean_OIC,resp_ORC,resp_OR_fractionC,resp_ANDC,resp_AND_fractionC] = get_responsive_fraction(RsC);

RsA = get_rasters_data(ei_A,selContexts,rasterNames);
mRsA = calc_mean_rasters(RsA,1:10);
RsA = find_responsive_rasters(RsA,1:10);
[CRcA,aCRcA,mRsRA] = find_population_vector_corr(RsA,mRsA,1);
% view_population_vector(Rs,mRs,400);
[resp_fractionA,resp_valsA,OIA,mean_OIA,resp_ORA,resp_OR_fractionA,resp_ANDA,resp_AND_fractionA] = get_responsive_fraction(RsA);
n = 0;

%%
%% population vector and correlation single animal
an = 4;
ff = makeFigureRowsCols(107,[1 0.5 4 1],'RowsCols',[2 4],...
    'spaceRowsCols',[0 -0.03],'rightUpShifts',[0.07 0.1],'widthHeightAdjustment',...
    [0.01 -60]);
set(gcf,'color','w');
set(gcf,'Position',[5 5 3.25 2]);
ff = show_population_vector_and_corr(mData,ff,RsC(an,:),mRsRC(an,:),CRcC(an,:));
save_pdf(ff.hf,mData.pdf_folder,sprintf('population_vector_corr.pdf'),600);

%% average correlation of all animals
ff = makeFigureRowsCols(107,[1 0.5 4 0.5],'RowsCols',[1 4],...
    'spaceRowsCols',[0 -0.03],'rightUpShifts',[0.075 0.2],'widthHeightAdjustment',...
    [0.01 -220]);
set(gcf,'color','w');
set(gcf,'Position',[5 5 3.25 1]);
ff = show_population_vector_and_corr(mData,ff,RsC(an,:),[],aCRcC);
save_pdf(ff.hf,mData.pdf_folder,sprintf('average_population_vector_corr.pdf'),600);

%% population vector and correlation single animal
an = 2;
ff = makeFigureRowsCols(107,[1 0.5 4 1],'RowsCols',[2 4],...
    'spaceRowsCols',[0 -0.03],'rightUpShifts',[0.07 0.1],'widthHeightAdjustment',...
    [0.01 -60]);
set(gcf,'color','w');
set(gcf,'Position',[5 5 3.25 2]);
ff = show_population_vector_and_corr(mData,ff,RsA(an,:),mRsRA(an,:),CRcA(an,:));
save_pdf(ff.hf,mData.pdf_folder,sprintf('population_vector_corr_AD.pdf'),600);

%% average correlation of all animals
ff = makeFigureRowsCols(107,[1 0.5 4 0.5],'RowsCols',[1 4],...
    'spaceRowsCols',[0 -0.03],'rightUpShifts',[0.075 0.2],'widthHeightAdjustment',...
    [0.01 -220]);
set(gcf,'color','w');
set(gcf,'Position',[5 5 3.25 1]);
ff = show_population_vector_and_corr(mData,ff,RsA(an,:),[],aCRcA);
save_pdf(ff.hf,mData.pdf_folder,sprintf('average_population_vector_corr_AD.pdf'),600);