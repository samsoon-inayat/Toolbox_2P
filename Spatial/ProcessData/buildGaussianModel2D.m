function [gmdl,x0] = buildGaussianModel2D(order,coeffs,dc)

inds = reshape(1:100,5,20)';
x0 = []; gmdl = '';
for ii = 1:order
    this_g = make_g_term2D(inds(ii,:));
    if ii == 1
        gmdl = this_g;
    else
        gmdl = [gmdl ' + ' this_g];
    end
    x0 = [x0,coeffs];
end
if dc
    dc_term = sprintf('b(%d)',inds(ii,5)+1);
    gmdl = [gmdl ' + ' dc_term];
    x0 = [x0,0];
end
gmdl = ['@(b,x)' gmdl];


function g_term = make_g_term2D(inds)
g_term = sprintf('b(%d)*exp(-((((x(:,1)-b(%d)).^2)/(2*(b(%d).^2))) + (((x(:,2)-b(%d)).^2)/(2*(b(%d).^2)))))',inds(1),inds(2),inds(3),inds(4),inds(5));
