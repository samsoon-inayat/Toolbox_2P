function view_population_vector(Rs,mRs,ccs,figNum)

if exist('figNum','var')
    figure(figNum);
else
    figure;
end

clf(gcf);
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
            if ~isempty(ccs)
                if iscell(ccs)
                    ccsr = ccs{rr,cc};
                else
                    if ccs == 1
                        ccsr = resp.vals;
                    else
                        ccsr = ccs;
                    end
                end
            else
                ccsr = logical(ones(size(1:length(resp.vals))));
            end
            [ptco,CRc,cellNums] = findPopulationVectorPlot(ptc,ccsr);
            imagesc(ptco);
            set(gca,'Ydir','Normal');
            colorbar;
            if rr == 1
                title(R.context_info);
            end
            if cc == 1
                ylabel(sprintf('Animal # %d',rr));
            end
            xlabel(sprintf('%.2f',sum(ccsr)/length(ccsr)));
        end
    end
else
end