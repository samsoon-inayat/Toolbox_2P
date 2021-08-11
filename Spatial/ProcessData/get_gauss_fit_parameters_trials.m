function [rs,as,bs,PWs] = get_gauss_fit_parameters_trials(mrfs,binwidth)

if isstruct(mrfs)
    % n = 0;
    rs = (squeeze(mrfs.coefficients_Rs_trials(:,4,:)))';
    A = (squeeze(mrfs.coefficients_Rs_trials(:,1,:)))';
    as = A;%*exp(0.5);
    bs = binwidth * (squeeze(mrfs.coefficients_Rs_trials(:,2,:)))';
    cs = (squeeze(mrfs.coefficients_Rs_trials(:,3,:)))';
    PWs = 2.36*cs./sqrt(2)*binwidth;
%     A = coeff(:,1);
%     mu = coeff(:,2);
%     sigma = coeff(:,3);
%     as = A';%*exp(0.5);
%     bs = mu'*binwidth;
%     cs = sigma';
%     PWs = 2.36*cs./sqrt(2)*binwidth;
else
%     for ii = 1:size(mrfs,1)
%             rs(ii) = mrfs(ii,4);
%             coeff(ii,:) = mrfs(ii,[1 2 3])';
%     end
%     A = coeff(:,1);
%     mu = coeff(:,2);
%     sigma = coeff(:,3);
%     as = A';%;*exp(0.5);
%     bs = mu'*binwidth;
%     cs = sigma';
%     PWs = 2.36*cs./sqrt(2)*binwidth;
end