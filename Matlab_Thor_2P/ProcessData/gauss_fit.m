function mSigF = gauss_fit(xs,coeffs,orderdc)
order = orderdc(1);
dc = orderdc(2);
[gmdl,~] = buildGaussianModel(order,coeffs,dc); mdl_fun_h = str2func(gmdl);
mSigF = mdl_fun_h(coeffs,xs);


% n = 0;