function view_population_vector_corr(Rs,mRs,figNum)

if exist('figNum','var')
    figure(figNum);
else
    figure;
end

if iscell(mRs)
    [nrows,ncols] = size(mRs);
    plotNums = 1:(nrows*ncols);
    plotNums = reshape(plotNums,ncols,nrows);
    plotNums = plotNums';
    for rr = 1:size(mRs,1)
        for cc = 1:size(mRs,2)
            subplot(nrows,ncols,plotNums(rr,cc));
            ptc = mRs{rr,cc};
            R = Rs{rr,cc};
            resp = R.resp;
            if ~isempty(resp)
                ccs = find(resp.vals);
            else
                ccs = [];
            end
            [ptco,CRc,cellNums] = findPopulationVectorPlot(ptc,ccs)
            imagesc(CRc);
            axis equal
            set(gca,'Ydir','Normal');
            colorbar;
            if rr == 1
                title(R.context_info);
            end
            if cc == 1
                ylabel(sprintf('Animal # %d',rr));
            end
            xlabel(sprintf('%.2f',resp.fraction));
        end
    end
else
end