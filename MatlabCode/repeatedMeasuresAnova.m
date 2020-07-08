function out = repeatedMeasuresAnova(data,varNames,within,varargin)

% p = inputParser;
% addRequired(p,'data',@isnumeric);
% addOptional(p,'dimension',1,@isnumeric);
% addOptional(p,'decimal_places',3,@isnumeric);
% addOptional(p,'do_mean','Yes');
% parse(p,data,varargin{:});
% 
% dimension = p.Results.dimension;
% decimal_places = p.Results.decimal_places;

cmdTxt = sprintf('dataT = table(');
for ii = 1:(size(data,2)-1)
    cmdTxt = sprintf('%sdata(:,%d),',cmdTxt,ii);
end
cmdTxt = sprintf('%sdata(:,size(data,2)));',cmdTxt);
eval(cmdTxt);
dataT.Properties.VariableNames = varNames;
withinVarNames = within.Properties.VariableNames;
for ii = 1:length(withinVarNames)
    cmdTxt = sprintf('within.%s = categorical(within.%s);',withinVarNames{ii},withinVarNames{ii});
    eval(cmdTxt);
end
withinModel = withinVarNames{1};
for ii = 2:length(withinVarNames)
    withinModel = [withinModel '*' withinVarNames{ii}];
end

% writetable(between,'Training_Data.xls');
cmdTxt = sprintf('rm = fitrm(dataT,''');
for ii = 1:(length(varNames)-1)
    cmdTxt = sprintf('%s%s,',cmdTxt,varNames{ii});
end
cmdTxt = sprintf('%s%s~1'');',cmdTxt,varNames{length(varNames)});
eval(cmdTxt);
rm.WithinDesign = within;
rm.WithinModel = withinModel;
out.rm = rm;
out.table = ranova(rm,'WithinModel',rm.WithinModel);
out.mauchlytbl = mauchly(rm);
% multcompare(rm,'Day','ComparisonType','bonferroni')
% mcTI = find_sig_mctbl(multcompare(rm,'TrialDiff','By','Condition','ComparisonType','bonferroni'),6);
% mcDays = find_sig_mctbl(multcompare(rm,'Condition','By','TrialDiff','ComparisonType','bonferroni'),6);
% [mVar semVar] = findMeanAndStandardError(data);
% 
% combs = nchoosek(1:size(data,2),2); p = ones(size(combs,1),1); h = logical(zeros(size(combs,1),1));
% for rr = 1:size(mcTI,1)
%     thisRow = mcTI(rr,:);
%     conditionN =  thisRow{1,1}; Rtype1 = thisRow{1,2}; Rtype2 = thisRow{1,3};
%     Num1 = find(ismember(within{:,:},[conditionN Rtype1],'rows'));
%     Num2 = find(ismember(within{:,:},[conditionN Rtype2],'rows'));
%     row = [Num1 Num2]; ii = ismember(combs,row,'rows'); p(ii) = mcTI{1,6}; h(ii) = 1;
% end
% for rr = 1:size(mcDays,1)
%     thisRow = mcDays(rr,:);
%     Rtype =  thisRow{1,1}; Condition1 = thisRow{1,2}; Condition2 = thisRow{1,3};
%     Num1 = find(ismember(within{:,:},[Condition1 Rtype],'rows'));
%     Num2 = find(ismember(within{:,:},[Condition2 Rtype],'rows'));
%     row = [Num1 Num2]; ii = ismember(combs,row,'rows'); p(ii) = mcDays{1,6}; h(ii) = 1;
% end