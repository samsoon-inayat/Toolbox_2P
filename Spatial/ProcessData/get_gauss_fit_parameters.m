function [rs,as,bs,PWs] = get_gauss_fit_parameters(mrfs,binwidth)
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

A = coeff(:,1);
mu = coeff(:,2);
sigma = coeff(:,3);
as = A'*exp(0.5);
bs = mu';
cs = sigma';
PWs = 2.36*cs./sqrt(2)*binwidth;