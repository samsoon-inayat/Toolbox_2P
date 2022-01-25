function [rs,as,bs,PWs] = get_gauss_fit_parameters(mrfs,binwidth)


if isstruct(mrfs)
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
    sigma = coeff(:,3); sigma = abs(sigma);
    as = A';%*exp(0.5);
    bs = mu'*binwidth;
    cs = sigma';
    PWs = 2.36*cs./sqrt(2)*binwidth;
else
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
end