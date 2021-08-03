function do_deconvolution(ops0,db)
ii = 1;
subFolder = [db(ii).mouse_name '\' db(ii).date '\' num2str(db(ii).expts)];
subFolderR = [db(ii).mouse_name '\' db(ii).date '\' num2str(db(ii).expText)];
pSaveDataFolder = makeName(subFolder,ops0.ResultsSavePath);
fileCheck = sprintf('%s\\F_*_plane1.mat',pSaveDataFolder);
files = dir(fileCheck);
fileName = makeName(files(1).name,pSaveDataFolder);
twoP = load(fileName);

caTraces = twoP.Fcell{1};

parfor ii = 1:length(twoP.stat)
    ii
%     if ii == 185
%         n = 0;
%     end
    meanCaTrace = mean(double(caTraces(ii,:)));
    caTrace = (double(caTraces(ii,:))'-meanCaTrace)/meanCaTrace;
    if any(isnan(caTrace))
        tsp = NaN(size(caTrace));
    else
    [ca_sig,tsp,options] = deconvolveCa(caTrace, 'ar1',...
    'constrained','optimize_b', 'optimize_pars', 'optimize_smin');
    if length(tsp) < length(caTrace)
        len_diff = length(caTrace) - length(tsp);
        tsp((length(tsp)+1):length(caTrace),1) = zeros(len_diff,1);
    end
    end
    sp(ii,:) = single(tsp');
end
twoP.p_sp = sp;
save(fileName,'-struct','twoP','-v7.3');


