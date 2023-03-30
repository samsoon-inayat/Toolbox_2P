function ppm_a = find_transients_per_minute(ei)

for ani = 1:length(ei)
    tei = ei{ani};
    if length(tei.plane) == 2
        tcas = [tei.plane{1}.tP.deconv.spSigAll;tei.plane{2}.tP.deconv.spSigAll];
    else
        tcas = tei.plane{1}.tP.deconv.spSigAll;
    end
    Fso = tei.thorExp.frameRate;
    T = 1/Fso;             % Sampling period       
    Lo = size(tcas,2);             % Length of signal
    t = (0:Lo-1)*T;        % Time vector
    minTime = min(t);
    maxTime = max(t);%maxTime = 15*60;
    ppm = NaN(size(tcas,1),1);
    parfor cni = 1:size(tcas,1)
        X = tcas(cni,:);
        [pks,locs,w,p] = findpeaks(X,Fso);
        ppm(cni) = length(pks)/(maxTime/60);
    end
    ppm_a{ani} = ppm;
end