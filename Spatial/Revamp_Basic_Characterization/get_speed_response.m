function [speedRs,resp_speed] = get_speed_response(ei)


var_names = {'linear','sigmoid','gauss'};
for ii = 1:length(ei)
    tei = ei{ii};
    psp = [];
    psp1 = tei.plane{1}.speed_response;
    if length(tei.plane) == 2
        psp2 = tei.plane{2}.speed_response;
        psp.corr = [psp1.corr;psp2.corr];
        psp.FR_vs_speed = [psp1.FR_vs_speed;psp2.FR_vs_speed];
         for vv = 1:length(var_names)
            cmdTxt = sprintf('psp.fits.%s.fitted = [psp1.fits.%s.fitted;psp2.fits.%s.fitted];',var_names{vv},var_names{vv},var_names{vv});
            eval(cmdTxt);
            cmdTxt = sprintf('psp.fits.%s.coeffsrs = [psp1.fits.%s.coeffsrs;psp2.fits.%s.coeffsrs];',var_names{vv},var_names{vv},var_names{vv});
            eval(cmdTxt);
         end
        psp.bin_centers = psp1.bin_centers;
        speedRs{ii,1} = psp;
    else
        speedRs{ii,1} = psp1;
    end
end
n = 0;
%% Percentage speed responsive

while 1
    for an = 1:5
        d.bcs = speedRs{an}.bin_centers;
        d.FR = speedRs{an}.FR_vs_speed;
        fitg = speedRs{an}.fits.gauss; fits = speedRs{an}.fits.sigmoid; fitl = speedRs{an}.fits.linear;
        d.fFRl = fitl.fitted; d.fFRs = fits.fitted; d.fFRg = fitg.fitted;
        d.cl = fitl.coeffsrs(:,3); d.cs = fits.coeffsrs(:,3); d.cg = fitg.coeffsrs(:,3);
        [rs,MFR,centers,PWs] = get_gauss_fit_parameters(fitg.coeffsrs,d.bcs(2)-d.bcs(1));
        inds = centers < 1 | centers > 39 | rs < 0.25 | PWs < 10;% | PWs > 20 | PWs < 10;
        inds = ~inds;
        pR(an) = 100*sum(inds)/length(inds);
        resp_speed{an,1} = inds';
    end
    break;
end