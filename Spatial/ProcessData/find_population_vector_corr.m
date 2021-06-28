function [CRc,aCRc,mR] = find_population_vector_corr(Rs,mRs,ccs,plt)

for rr = 1:size(mRs,1)
    for cc = 1:size(mRs,2)
        ptc = mRs{rr,cc};
        R = Rs{rr,cc};
        resp = R.resp;
        if ~isempty(ccs)
            if ccs == 1
                ccsr = find(resp.vals);
            end
        else
            ccsr = ccs;
        end
        [ptco,temp_C,~] = findPopulationVectorPlot(ptc,ccsr);
%         [~,temp_C,~] = findPopulationVectorPlot(ptc,ccsr);
        CRc{rr,cc} = temp_C;
        sz_CRc(rr,cc) = size(temp_C,1);
        mR{rr,cc} = ptco;
   end
end

for cc = 1:size(mRs,2)
    msz = max(sz_CRc(:,cc));
    temp = NaN(msz,msz);
    aCRc{cc} = repmat(temp,1,1,size(mRs,1));
end


for cc = 1:size(mRs,2)
    for rr = 1:size(mRs,1)
        temp_C = CRc{rr,cc};
        aCRc{cc}(1:size(temp_C,1),1:size(temp_C,2),rr) = temp_C;
    end
    aCRc{cc} = nanmean(aCRc{cc},3);
end

if exist('plt','var')
    figure(plt);clf
    for cc = 1:size(mRs,2)
        subplot(1,size(mRs,2),cc);
        imagesc(aCRc{cc});
        axis equal
%         axis off;
        set(gca,'Ydir','Normal');
%         colorbar;
    end
end

