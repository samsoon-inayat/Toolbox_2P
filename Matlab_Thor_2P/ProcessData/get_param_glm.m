function [values,cns,acs] = get_param_glm(paramName,ei,glms,an,contextN)

values = [];
cns = [];
acs = [];
nplanes = length(ei{an}.plane);
for pp = 1:nplanes
    aglm = glms{an,pp,contextN};
    temp = aglm;%get_val_glm(aglm);
    values = combineData(values,temp);
    cns = [cns;ones(size(temp.STD_S,2),1)*an,ones(size(temp.STD_S,2),1)*pp,(1:size(temp.STD_S,2))',ei{an}.plane{pp}.tP.iscell(:,1)];
end


function out = get_val_glm(aglm)

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

for cn = 1:size(aglm.mdl1,1)
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
        cmdTxt = sprintf('out.%s(cn) = 2*(aglm.mdl%d{cn}.LogLikelihood - aglm.mdl%d{cn}.LogLikelihood);',this,m1,m2);
        eval(cmdTxt);
    end
end
n = 0;

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
