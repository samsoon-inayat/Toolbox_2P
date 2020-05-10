function [rs,coeff] = getMRFS_vals(mrfs)
% n = 0;
for ii = 1:length(mrfs.worked)
    if mrfs.worked(ii)
        rs(ii) = mrfs.coefficients_Rs_mean(ii,4);
        coeff(ii,:) = mrfs.coefficients_Rs_mean(ii,[1 2 3])';
    else
        rs(ii) = NaN;
        coeff(ii,:) = [NaN NaN NaN];
    end
end