function [CRc,aCRc,mRR,spCR] = find_population_correlations(Rs,mR,resp,ordercol)

[CRc,aCRc,mRR] = find_population_vector_corr(Rs,mR,resp,ordercol);

for rr = 1:size(mRR,1)
    for cc = 1:size(mRR,2)
        thismR = mRR{rr,cc};
%         thismR = thismR(Rs{rr,cc}.resp.vals,:);
        spCR{rr,cc} = corr(thismR');
    end
end

n = 0;