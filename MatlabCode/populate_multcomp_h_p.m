function [combs,h,p] = populate_multcomp_h_p(all_data,within,mcTI,mcDays)

allRows = categorical(within{:,:});

combs = nchoosek(1:size(all_data,2),2); p = ones(size(combs,1),1); h = logical(zeros(size(combs,1),1));
for rr = 1:size(mcTI,1)
    thisRow = mcTI(rr,:);
    conditionN =  thisRow{1,1}; Rtype1 = thisRow{1,2}; Rtype2 = thisRow{1,3};
    Num1 = find(ismember(allRows,[conditionN Rtype1],'rows'));
    Num2 = find(ismember(allRows,[conditionN Rtype2],'rows'));
    row = [Num1 Num2]; ii = ismember(combs,row,'rows'); p(ii) = mcTI{1,6}; h(ii) = 1;
end
for rr = 1:size(mcDays,1)
    thisRow = mcDays(rr,:);
    Rtype =  thisRow{1,1}; Condition1 = thisRow{1,2}; Condition2 = thisRow{1,3};
    Num1 = find(ismember(allRows,[Condition1 Rtype],'rows'));
    Num2 = find(ismember(allRows,[Condition2 Rtype],'rows'));
    row = [Num1 Num2]; ii = ismember(combs,row,'rows'); p(ii) = mcDays{1,6}; h(ii) = 1;
end