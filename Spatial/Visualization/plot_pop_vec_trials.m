function [ff,seq] = plot_pop_vec_trials(FN,pos,tRs,tmR,resp,seqi)

mData = evalin('base','mData');

respT = repmat(resp,1,size(tmR,2));

ff = makeFigureRowsCols(FN,[1 0.5 4 0.5],'RowsCols',[2 11],...
        'spaceRowsCols',[0.02 0.025],'rightUpShifts',[0.03 0.13],'widthHeightAdjustment',...
        [-30 -200]);    set(gcf,'color','w');    set(gcf,'Position',pos);
    if ~exist('seqi','var')
        [CRc,aCRc,mRR,seq] = find_population_vector_corr(tRs,tmR,respT,1);
    else
        [CRc,aCRc,mRR,seq] = find_population_vector_corr(tRs,tmR,respT,1,seqi);
    end
%     ff = show_population_vector_and_corr_peak(mData,ff,tRs,mRR,CRc,[],[]);
    ff = show_population_vector_and_corr(mData,ff,tRs,mRR,CRc,[],[]);