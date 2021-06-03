function [CRc,aCRc] = find_population_vector_corr(Rs,mRs,plt)

for rr = 1:size(mRs,1)
    for cc = 1:size(mRs,2)
        ptc = mRs{rr,cc};
        R = Rs{rr,cc};
        resp = R.resp;
        if ~isempty(resp)
            ccs = find(resp.vals);
        else
            ccs = [];
        end
        [~,CRc{rr,cc},~] = findPopulationVectorPlot(ptc,ccs)
   end
end

for cc = 1:size(mRs,2)
    temp = zeros(size(CRc{1,cc}));
    aCRc{cc} = repmat(temp,1,1,size(mRs,1));
end

for cc = 1:size(mRs,2)
    for rr = 1:size(mRs,1)
        aCRc{cc}(:,:,rr) = CRc{rr,cc};
    end
    aCRc{cc} = nanmean(aCRc{cc},3);
end

if plt(1) == 1
    figure(500);clf
    for cc = 1:size(mRs,2)
        subplot(1,size(mRs,2),cc);
        imagesc(aCRc{cc});
        axis equal
        axis off;
        set(gca,'Ydir','Normal');
%         colorbar;
    end
end

