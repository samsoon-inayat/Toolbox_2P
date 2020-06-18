function [rs,coeff] = getMRFS_vals_Trials(mrfs)
% n = 0;
coeff = [];
for ii = 1:length(mrfs.worked)
    rs(:,ii) = mrfs.coefficients_Rs_trials(:,4,ii);
    coeff(:,:,ii) = mrfs.coefficients_Rs_trials(:,1:3,ii);
end
