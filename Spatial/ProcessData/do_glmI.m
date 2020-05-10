function glms = do_glm(aei,owr)

if ~exist('aei','var')
    aei = evalin('base','ei10(1)');
    owr = 1;
end
for ii = 1:length(aei)
    ei = aei{ii};
    for pp = 1:length(ei.plane)
        tplane = ei.plane{pp};
        for cc = 1:length(tplane.contexts)
            disp(sprintf('%s - Plane %d - Context %d',tplane.folder,pp,cc));
            fileName = fullfile(tplane.folder,sprintf('glmsI_%d.mat',cc));
            if ~exist(fileName,'file') || owr == 1
                glms{ii,pp,cc} = do_glm1(ei,pp,cc,fileName);
            else
                glms{ii,pp,cc} = load(fileName);
            end
        end
    end
end

function out = do_glm1(ei,planeNumber,contextNumber,fileName)
% planeNumber = 1;
% contextNumber = 1;
maxDistTime = [Inf Inf];

[dataTi,~,~] = getParamValues('',{ei},planeNumber,contextNumber,'airI','time','All',maxDistTime);
[dataD,~,~] = getParamValues('',{ei},planeNumber,contextNumber,'air','dist','All',maxDistTime);
dataT = dataTi.fromFrames; tempSampl = diff(dataT.duration(1,:)); samplingRate = tempSampl(1); allTs = (0:samplingRate:50);
dataT.xs = allTs(1:size(dataT.sp_rasters,2));
modelName = {'STD';'ST';'TD';'SD';'S';'T';'D'};
devianceVars = {'STD_S';'STD_TD';'STD_T';'STD_D';'ST_T';'ST_S';'SD_D';'SD_S'};
modelTermsS = {'x'};
modelTermsT = {'t + t^2 + t^3 + t^4 + t^5'};
modelTermsD = {'d + d^2 + d^3 + d^4 + d^5'};
% modelTermsCommon = {'FR ~ 1 + s + h1 + h2 + h3 + h4 + h5 + h6'};
modelTermsCommon = {'FR ~ 1 + s + h5 + h6'};
modelspec = cell(length(modelName),1);
for ii = 1:length(modelName)
    modelspec{ii} = modelTermsCommon;
    this = modelName{ii};
    if strfind(this,'S')
        modelspec{ii} = strcat(modelspec{ii},' + ',modelTermsS);
    end
    if strfind(this,'T')
        modelspec{ii} = strcat(modelspec{ii},' + ',modelTermsT);
    end
    if strfind(this,'D')
        modelspec{ii} = strcat(modelspec{ii},' + ',modelTermsD);
    end
end
% dataT.sp_rasters = dataT.sp_rasters_nan_corrected;
firing_history = construct_firing_history(dataT);
timeVar = repmat(dataT.xs(1:size(dataT.sp_rasters,2))',size(dataT.sp_rasters,1),1);
distVar = reshape(dataT.dist(:,1:size(dataT.sp_rasters,2))',numel(dataT.sp_rasters(:,:,1)),1);
spaceVar = reshape(dataT.space(:,1:size(dataT.sp_rasters,2))',numel(dataT.sp_rasters(:,:,1)),1);
speedVar = reshape(dataT.speed(:,1:size(dataT.sp_rasters,2))',numel(dataT.sp_rasters(:,:,1)),1);

mdl1 = cell(size(firing_history,3),1);
mdl2 = mdl1; mdl3 = mdl1; mdl4 = mdl1; mdl5 = mdl1; mdl6 = mdl1; mdl7 = mdl1;
mdl1s = mdl1; mdl2s = mdl1; mdl3s = mdl1; mdl4s = mdl1; mdl5s = mdl1; mdl6s = mdl1; mdl7s = mdl1;
parfor cn = 1:size(firing_history,3)
%     h1 = firing_history(:,1,cn);     h2 = firing_history(:,2,cn);     h3 = firing_history(:,3,cn); 
%     h4 = firing_history(:,4,cn);     
    h5 = firing_history(:,5,cn);     h6 = firing_history(:,6,cn); 
    ResponseVar = reshape(dataT.sp_rasters(:,:,cn),numel(dataT.sp_rasters(:,:,cn)),1);
%     data_table = table(timeVar,distVar,spaceVar,speedVar,h1,h2,h3,h4,h5,h6,ResponseVar,'VariableNames',{'t','d','x','s','h1','h2','h3','h4','h5','h6','FR'});
    data_table = table(timeVar,distVar,spaceVar,speedVar,h5,h6,ResponseVar,'VariableNames',{'t','d','x','s','h5','h6','FR'});
    mdl1{cn} = fitglm(data_table,modelspec{1}{1},'Distribution','Poisson'); %mdl = step(mdl1{1},'NSteps',50);
    mdl2{cn} = fitglm(data_table,modelspec{2}{1},'Distribution','Poisson');
    mdl3{cn} = fitglm(data_table,modelspec{3}{1},'Distribution','Poisson');
    mdl4{cn} = fitglm(data_table,modelspec{4}{1},'Distribution','Poisson');
    mdl5{cn} = fitglm(data_table,modelspec{5}{1},'Distribution','Poisson');
    mdl6{cn} = fitglm(data_table,modelspec{6}{1},'Distribution','Poisson');
    mdl7{cn} = fitglm(data_table,modelspec{7}{1},'Distribution','Poisson');
end
n = 0;

for cn = 1:size(firing_history,3)
    temp = mdl1{cn}.devianceTest;
    val = temp;
end

for cn = 1:size(firing_history,3)
    for ii = 1:length(devianceVars)
        this = devianceVars{ii};
        ind = strfind(this,'_');
        first = this(1:(ind-1)); second = this((ind+1):end);
        ind = strcmp(modelName,first);
        m1 = find(ind);
        diffstr = setdiff(first,second);
        ind = strcmp(modelName,diffstr);
        if sum(ind) == 0
            ind = strcmp(modelName,fliplr(diffstr));
        end
        m2 = find(ind);
%         disp([this ' ' modelName{m1} ' ' modelName{m2}])
        cmdTxt = sprintf('out.%s(cn) = 2*(mdl%d{cn}.LogLikelihood - mdl%d{cn}.LogLikelihood);',this,m1,m2);
        eval(cmdTxt);
        nterms1 = length(strfind(modelspec{m1}{1},'+'))+1;
        nterms2 = length(strfind(modelspec{m2}{1},'+'))+1;
        cmdTxt = sprintf('out.pvals(cn,ii) = 1 - chi2cdf(out.%s(cn),nterms1-nterms2);',this);
        eval(cmdTxt);
    end
end
n = 0;

% fileName = fullfile(ei.plane{planeNumber}.folder,sprintf('glms_%d.mat',contextNumber));
% for ii = 1:length(modelspec)
%     cmdTxt = sprintf('out.mdl%d = mdl%d;',ii,ii);
%     eval(cmdTxt);
% end
save(fileName,'-struct','out','-v7.3');


% 
% d_Dev = STD_T - STD_D;
% d_Dev = SD_S - SD_D;
% figure(10000);clf;hist(d_Dev,50);
% 
% figure(10000);clf;
% scatter(STD_T,STD_D);
% set(gca,'XScale','log','YScale','log');
% 
% 
% cns = find(d_Dev > 100);
% for cc = 1:length(cns)
%     cn = cns(cc);
%     figure(2000);clf;
%     subplot 121;
%     imagesc(dataT.sp_rasters(:,:,cn));title('Time');
%     subplot 122;
%     imagesc(dataD.sp_rasters(:,:,cn));title('Dist');
%     pause;
% end




function firing_history = construct_firing_history(dataT)
firing_history = NaN(numel(dataT.sp_rasters(:,:,1)),size(dataT.cell_history,2));
for cn = 1:size(dataT.cell_history,3)
    thisRaster = dataT.sp_rasters(:,:,cn);
    cellH = dataT.cell_history(:,:,cn);
    bigRaster = [cellH thisRaster];
    for nn = 1:6
        thisFiringHistory = bigRaster(:,nn:(nn+size(thisRaster,2)-1));
        firing_history(:,nn,cn) = reshape(thisFiringHistory,numel(thisFiringHistory),1);
    end
end