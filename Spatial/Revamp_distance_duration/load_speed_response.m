function [speedRs,resp_speed,pR,resp_speedN] = get_speed_response(ei)


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
        d.bcs = speedRs{an,1}.bin_centers;
        d.FR = speedRs{an,1}.FR_vs_speed;
        fitg = speedRs{an,1}.fits.gauss; fits = speedRs{an,1}.fits.sigmoid; fitl = speedRs{an,1}.fits.linear;
        d.fFRl = fitl.fitted; d.fFRs = fits.fitted; d.fFRg = fitg.fitted;
        d.cl = fitl.coeffsrs(:,3); d.cs = fits.coeffsrs(:,3); d.cg = fitg.coeffsrs(:,3);
        [rs,MFR,centers,PWs] = get_gauss_fit_parameters(fitg.coeffsrs,d.bcs(2)-d.bcs(1));
        inds = centers < 1 | centers > 39 | rs < 0.3 | PWs < 10;% | PWs > 20 | PWs < 10;
        inds = centers < 1 | centers > 39 | PWs < 1 | PWs > 40 | rs < 0.3;
        inds = ~inds;
        pR(an) = 100*sum(inds)/length(inds);
        resp_speed{an,1} = inds';
    end
    break;
end
%%
for ii = 1:length(ei)
    tei = ei{ii};
    psp = [];
    psp1 = tei.plane{1}.speed_response.McN;
    if length(tei.plane) == 2
        psp2 = tei.plane{2}.speed_response.McN;
        psp.speed_resp = [psp1.speed_resp;psp2.speed_resp];
        psp.speed_tuning_inc = [psp1.speed_tuning_inc;psp2.speed_tuning_inc];
        psp.speed_tuning_dec = [psp1.speed_tuning_dec;psp2.speed_tuning_dec];
        
        psp.bins = psp1.bins;
        speedRs{ii,2} = psp;
        resp_speed{ii,2} = psp.speed_resp;
        resp_speed{ii,3} = abs(psp.speed_resp);
    else
        speedRs{ii,2} = psp1;
        resp_speed{ii,2} = psp1.speed_resp;
        resp_speed{ii,3} = logical(abs(psp1.speed_resp));
    end
end

%%
var_names = {'gauss'};
for ii = 1:length(ei)
    tei = ei{ii};
    psp = [];
    psp1 = tei.plane{1}.speed_response_gauss;
    if length(tei.plane) == 2
        psp2 = tei.plane{2}.speed_response_gauss;
        psp.corr = [psp1.corr;psp2.corr];
        psp.corra = [psp1.corra;psp2.corra];
        psp.FR_vs_speed = [psp1.FR_vs_speed;psp2.FR_vs_speed];
        psp.FR_vs_accel = [psp1.FR_vs_accel;psp2.FR_vs_accel];
        psp.fits.gaussR.fitted = [psp1.fits.gaussR.fitted;psp2.fits.gaussR.fitted];
        psp.fits.gaussR.coeffsrs = [psp1.fits.gaussR.coeffsrs;psp2.fits.gaussR.coeffsrs];
        psp.fits.gaussR.fittedS = cat(1,psp1.fits.gaussR.fittedS,psp2.fits.gaussR.fittedS);
        psp.fits.gaussR.coeffsrsS= cat(1,psp1.fits.gaussR.coeffsrsS,psp2.fits.gaussR.coeffsrsS);

        psp.fitsA.gaussR.fitted = [psp1.fitsA.gaussR.fitted;psp2.fitsA.gaussR.fitted];
        psp.fitsA.gaussR.coeffsrs = [psp1.fitsA.gaussR.coeffsrs;psp2.fitsA.gaussR.coeffsrs];
        psp.fitsA.gaussR.fittedS = cat(1,psp1.fitsA.gaussR.fittedS,psp2.fitsA.gaussR.fittedS);
        psp.fitsA.gaussR.coeffsrsS= cat(1,psp1.fitsA.gaussR.coeffsrsS,psp2.fitsA.gaussR.coeffsrsS);

        psp.bin_centers = psp1.bin_centers;
        psp.bin_centersA = psp1.bin_centersA;
        speedRs{ii,3} = psp;
    else
        speedRs{ii,3} = psp1;
    end
end
n = 0;

%% speed gauss find responsiveness
while 1
    speedRG = speedRs(:,3);
    for ii = 1:length(speedRG)
        tsrg = speedRG{ii};
        gaussR = tsrg.fits.gaussR;
        bcs = tsrg.bin_centers;
        [rs,MFR,centers,PWs] = get_gauss_fit_parameters(gaussR.coeffsrs,bcs(2)-bcs(1));
        nshuffle = size(gaussR.coeffsrsS,3);
        clear rsS;
        for jj = 1:nshuffle
            tcoeff = gaussR.coeffsrsS(:,:,jj);
            [rsS(jj,:),~,~,~] = get_gauss_fit_parameters(tcoeff,bcs(2)-bcs(1));
        end
        rsST = rsS';
        rs_r = repmat(rs',1,nshuffle);
        p_vals = sum(rs_r > rsST,2)/nshuffle;
        sp_cells = p_vals > 0.95;
%         for cn = 1:length(sp_cells)
%             if sp_cells(cn) == 0
%                 continue;
%             end
%             figure(10000);clf;
%             plot(bcs,tsrg.FR_vs_speed(cn,:));hold on;
%             plot(bcs,gaussR.fitted(cn,:));
%             pause(0.1);
%         end
        inds = centers < bcs(1) | centers > bcs(end) | PWs > (bcs(end)-bcs(1));
%         inds = rs < 0.7;
        resp_speed{ii,4} = sp_cells & ~(inds');
        resp_speedN{ii,1} = resp_speed{ii,4};
    end
    %%
    break;
end

%% accel gauss find responsiveness
while 1
    speedRG = speedRs(:,3);
    for ii = 1:length(speedRG)
        tsrg = speedRG{ii};
        gaussR = tsrg.fitsA.gaussR;
        bcs = tsrg.bin_centersA;
        [rs,MFR,centers,PWs] = get_gauss_fit_parameters(gaussR.coeffsrs,bcs(2)-bcs(1));
        nshuffle = size(gaussR.coeffsrsS,3);
        clear rsS;
        for jj = 1:nshuffle
            tcoeff = gaussR.coeffsrsS(:,:,jj);
            [rsS(jj,:),~,~,~] = get_gauss_fit_parameters(tcoeff,bcs(2)-bcs(1));
        end
        rsST = rsS';
        rs_r = repmat(rs',1,nshuffle);
        p_vals = sum(rs_r > rsST,2)/nshuffle;
        sp_cells = p_vals > 0.95;
%         for cn = 1:length(sp_cells)
%             if PWs(cn)>0
%                 continue;
%             end
%             figure(10000);clf;
%             plot(bcs,tsrg.FR_vs_speed(cn,:));hold on;
%             plot(bcs,gaussR.fitted(cn,:));
%             pause(0.1);
%         end
        inds = centers < bcs(1) | centers > bcs(end) | PWs > (bcs(end)-bcs(1));
%         inds = rs < 0.7;
        resp_speed{ii,5} = sp_cells & ~(inds');
        resp_speedN{ii,2} = resp_speed{ii,5};
    end
    %%
    break;
end


var_names = {'linear','sigmoid','gauss'};
for ii = 1:length(ei)
    tei = ei{ii};
    psp = [];
    psp1 = tei.plane{1}.accel_response;
    if length(tei.plane) == 2
        psp2 = tei.plane{2}.accel_response;
        psp.corra = [psp1.corra;psp2.corra];
        psp.FR_vs_accel = [psp1.FR_vs_accel;psp2.FR_vs_accel];
         for vv = 1:length(var_names)
            cmdTxt = sprintf('psp.fits.%s.fitted = [psp1.fits.%s.fitted;psp2.fits.%s.fitted];',var_names{vv},var_names{vv},var_names{vv});
            eval(cmdTxt);
            cmdTxt = sprintf('psp.fits.%s.coeffsrs = [psp1.fits.%s.coeffsrs;psp2.fits.%s.coeffsrs];',var_names{vv},var_names{vv},var_names{vv});
            eval(cmdTxt);
         end
        psp.bin_centers = psp1.bin_centers;
        speedRs{ii,4} = psp;
    else
        speedRs{ii,4} = psp1;
    end
end
n = 0;
%%  accel responsive
while 1
    for an = 1:5
        d.bcs = speedRs{an,4}.bin_centers;
        d.FR = speedRs{an,4}.FR_vs_accel;
        fitg = speedRs{an,4}.fits.gauss; fits = speedRs{an,4}.fits.sigmoid; fitl = speedRs{an,4}.fits.linear;
        d.fFRl = fitl.fitted; d.fFRs = fits.fitted; d.fFRg = fitg.fitted;
        d.cl = fitl.coeffsrs(:,3); d.cs = fits.coeffsrs(:,3); d.cg = fitg.coeffsrs(:,3);
        [rs,MFR,centers,PWs] = get_gauss_fit_parameters(fitg.coeffsrs,d.bcs(2)-d.bcs(1));
        inds = centers < 1 | centers > 39 | rs < 0.3 | PWs < 10;% | PWs > 20 | PWs < 10;
        inds = centers < 1 | centers > 39 | PWs < 1 | PWs > 40 | rs < 0.3;
        inds = ~inds;
        pR(an) = 100*sum(inds)/length(inds);
        resp_speed{an,6} = inds';
    end
    break;
end