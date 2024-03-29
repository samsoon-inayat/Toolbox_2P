function [resp_fraction,resp_vals,all_OI,mean_OI,resp_OR,resp_OR_fraction,resp_AND,resp_AND_fraction,resp_exc_inh] = get_responsive_fraction(Rs)

% resp_exc_inh = NaN;

for rr = 1:size(Rs,1)
    ccs = [];
    exc = [];
    inh = [];
    for cc = 1:size(Rs,2)
        R = Rs{rr,cc};
        try
            ccs(:,cc) = R.resp.vals';
        catch
            ccs(:,cc) = R.resp.valsA';
        end
        try
            resp_fraction(rr,cc) = R.resp.fraction;
        catch
            resp_fraction(rr,cc) = R.resp.fractionA;
        end
        if isfield(R.resp,'excinh')
            exc(:,cc) = R.resp.excinh==1;
            inh(:,cc) = R.resp.excinh==0;
        else
            exc_inh = NaN;
        end
    end
    resp_vals{rr} = logical(ccs);
    resp_exc_inh{rr,1} = logical(exc);
    resp_exc_inh{rr,2} = logical(inh);
    
    resp_OR{rr} = logical(zeros(size(ccs,1),1));
    resp_AND{rr} = logical(ones(size(ccs,1),1));
    for cc = 1:size(Rs,2)
        resp_OR{rr} = resp_OR{rr} | logical(ccs(:,cc));
        resp_AND{rr} = resp_AND{rr} & logical(ccs(:,cc));
    end
    resp_OR_fraction(rr) = sum(resp_OR{rr})/length(resp_OR{rr});
    resp_AND_fraction(rr) = sum(resp_AND{rr})/length(resp_AND{rr});
end

for ii = 1:length(resp_vals)
    ccs = resp_vals{ii};
    OI = NaN(size(ccs,2),size(ccs,2));
    mask = triu(ones(size(ccs,2)),1);
    for rr = 1:size(ccs,2)
        ccs1 = ccs(:,rr);
        for cc = 1:size(ccs,2)
            if mask(rr,cc)
                ccs2 = ccs(:,cc);
                shared = ccs1 & ccs2;
                OI(rr,cc) = sum(shared)/(sum(ccs1)+sum(ccs2)-sum(shared));
            end
        end
    end
    all_OI{ii} = OI;
end

all_OIm = repmat(OI,1,1,size(ccs,2));
for ii = 1:length(resp_vals)
    OI = all_OI{ii};
    all_OIm(:,:,ii) = OI;
end

mean_OI = nanmean(all_OIm,3);