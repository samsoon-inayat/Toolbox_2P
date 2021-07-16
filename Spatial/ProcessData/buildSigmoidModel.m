function [gmdl,x0] = buildSigmoidModel(order,coeffs,dc)

inds = reshape(1:98,2,49)';
x0 = []; gmdl = '';
for ii = 1:order
    this_g = make_g_term(inds(ii,:));
    if ii == 1
        gmdl = this_g;
    else
        gmdl = [gmdl ' + ' this_g];
    end
    x0 = [x0,coeffs];
end
if dc
    dc_term = sprintf('b(%d)',inds(ii,3)+1);
    gmdl = [gmdl ' + ' dc_term];
    x0 = [x0,0];
end
gmdl = ['@(b,x)' gmdl];


function g_term = make_g_term(inds)
g_term = sprintf('1./(1+exp(-b(%d)*(x-b(%d))))',inds(1),inds(2));