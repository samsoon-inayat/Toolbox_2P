function [rs,as,bs,PWs] = get_gauss_fit_parameters(mrfs,bins)


for ii = 1:size(mrfs,1)
        rs(ii) = mrfs(ii,4);
        coeff(ii,:) = mrfs(ii,[1 2 3])';
end
A = coeff(:,1);
mu = coeff(:,2);
sigma = coeff(:,3); sigma = abs(sigma);
as = A';%;*exp(0.5);
bs = mu'*binwidth;
cs = sigma';
PWs = 2.36*cs./sqrt(2)*binwidth;
