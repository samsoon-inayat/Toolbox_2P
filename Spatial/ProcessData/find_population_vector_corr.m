function [CRc,aCRc,mR,seq,msz] = find_population_vector_corr(Rs,mRs,ccs,ordercol,seqi)


for rr = 1:size(mRs,1)
    for cc = 1:size(mRs,2)
        R = Rs{rr,cc};
        resp = R.resp;
        if ~iscell(ccs)
            if ccs == 0
                accsr{rr,cc} = [];
            else
                accsr{rr,cc} = find(resp.vals);
            end
        else
            accsr{rr,cc} = find(ccs{rr,cc});
        end
        ccsr = accsr{rr,cc};
        ptc = mRs{rr,cc};
        [~,~,cellNums{rr,cc}] = findPopulationVectorPlot(ptc,ccsr);
    end
end

for rr = 1:size(mRs,1)
    for cc = 1:size(mRs,2)
        ptc = mRs{rr,cc};
        R = Rs{rr,cc};
        resp = R.resp;
        ccsr = accsr{rr,cc};
        if ~exist('seqi','var')
        if ordercol > 0
            try
                [ptco,temp_C,~] = findPopulationVectorPlot(ptc,ccsr,cellNums{rr,ordercol});
            catch
                disp('To order properly the number of cells must be same across conditions');
                lasterror
            end
        else
            [ptco,temp_C,~] = findPopulationVectorPlot(ptc,ccsr);
        end
        else
            [ptco,temp_C,~] = findPopulationVectorPlot(ptc,ccsr,seqi);
        end
        CRc{rr,cc} = temp_C;
        sz_CRc(rr,cc) = size(temp_C,1);
        mR{rr,cc} = ptco;
   end
end

for cc = 1:size(mRs,2)
    msz(cc) = max(sz_CRc(:,cc));
    temp = NaN(msz(cc),msz(cc));
    aCRc{cc} = repmat(temp,1,1,size(mRs,1));
end


for cc = 1:size(mRs,2)
    for rr = 1:size(mRs,1)
        temp_C = CRc{rr,cc};
        aCRc{cc}(1:size(temp_C,1),1:size(temp_C,2),rr) = temp_C;
    end
    aCRc{cc} = nanmean(aCRc{cc},3);
end
seq = cellNums;
n = 0;
% if exist('plt','var')
%     figure(plt);clf
%     for cc = 1:size(mRs,2)
%         subplot(1,size(mRs,2),cc);
%         imagesc(aCRc{cc});
%         axis equal
% %         axis off;
%         set(gca,'Ydir','Normal');
% %         colorbar;
%     end
% end

