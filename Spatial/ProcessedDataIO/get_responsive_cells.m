function out = get_responsive_cells(Rs)
% function [resp_fraction,resp_vals,all_OI,mean_OI,resp_OR,resp_OR_fraction,resp_AND,resp_AND_fraction,resp_exc_inh] = get_responsive_cells(Rs)
for rr = 1:size(Rs,1)
    for cc = 1:size(Rs,2)
        R = Rs{rr,cc};
        out.vals{rr,cc} = R.resp.vals';
        out.fraction(rr,cc) = R.resp.fraction;
        if isfield(R.resp,'excinh')
            out.exc(:,cc) = R.resp.excinh==1;
            out.inh(:,cc) = R.resp.excinh==0;
        else
            out.exc(:,cc) = NaN;
            out.inh(:,cc) = NaN;
        end
    end
end
