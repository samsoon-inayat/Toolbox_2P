function do_glm(ei)

if ~exist('ei','var')
    aei = evalin('base','ei10');
    ei = aei{7};
end



planeNumber = 1;
contextNumber = 3;
maxDistTime = [142 5];

[dataT1 cns areCells] = getParamValues('',{ei},planeNumber,contextNumber,'air','time','All',maxDistTime);
dataT = dataT1.wholeData;
dataT.sp_rasters_nan_corrected = dataT.sp_rasters;
firing_history = construct_firing_history(dataT);
[dataD cns areCells] = getParamValues('',{ei},planeNumber,contextNumber,'air','dist','All',maxDistTime);


modelName = {'STD';'ST';'TD';'SD';'S';'T';'D'};
devianceVars = {'STD_S';'STD_TD';'STD_T';'STD_D';'ST_T';'ST_S';'SD_D';'SD_S'};
modelTermsS = {'x'};
modelTermsT = {'t + t^2 + t^3 + t^4 + t^5'};
modelTermsD = {'d + d^2 + d^3 + d^4 + d^5'};
% modelTermsCommon = {'FR ~ 1 + s + h1 + h2 + h3 + h4 + h5 + h6'};
modelTermsCommon = {'FR ~ 1 + s + h5 + h6'};
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


timeVar = dataT.xs';
distVar = dataT.dist;
spaceVar = dataT.space;
speedVar = dataT.speed;

for cn = 1:size(firing_history,2)
    cn
    for ii = 1:size(firing_history,3)
        cmdTxt = sprintf('h%d = firing_history(:,cn,ii);',ii);
        eval(cmdTxt);
    end
    ResponseVar = dataT.sp_rasters_nan_corrected(:,cn);
%     data_table = table(timeVar,distVar,spaceVar,speedVar,h1,h2,h3,h4,h5,h6,ResponseVar,'VariableNames',{'t','d','x','s','h1','h2','h3','h4','h5','h6','FR'});
    data_table = table(timeVar,distVar,spaceVar,speedVar,h5,h6,ResponseVar,'VariableNames',{'t','d','x','s','h5','h6','FR'});
    for ii = 1:length(modelName)
        mdl{ii} = fitglm(data_table,modelspec{ii}{1},'Distribution','Poisson');
    end
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
        cmdTxt = sprintf('%s(cn) = 2*(mdl{%d}.LogLikelihood - mdl{%d}.LogLikelihood);',this,m1,m2);
        eval(cmdTxt);
    end
end
n = 0;

storeVars = {'modelName';'devianceVars'

d_Dev = STD_T - STD_D;
d_Dev = SD_S - SD_D;
figure(10000);clf;hist(d_Dev,50);

figure(10000);clf;
scatter(STD_T,STD_D);
set(gca,'XScale','log','YScale','log');


cns = find(d_Dev > 100);
for cc = 1:length(cns)
    cn = cns(cc);
    figure(2000);clf;
    subplot 121;
    imagesc(dataT1.sp_rasters(:,:,cn));title('Time');
    subplot 122;
    imagesc(dataD.sp_rasters(:,:,cn));title('Dist');
    pause;
end




function firing_history = construct_firing_history(dataT)
firing_history = NaN(size(dataT.sp_rasters_nan_corrected,1),size(dataT.sp_rasters_nan_corrected,2),size(dataT.cell_history,2));
bigRaster = [dataT.cell_history';dataT.sp_rasters_nan_corrected];
for nn = 1:6
    thisFiringHistory = bigRaster(nn:(nn+size(dataT.sp_rasters_nan_corrected,1)-1),:);
    firing_history(:,:,nn) = thisFiringHistory;
end
