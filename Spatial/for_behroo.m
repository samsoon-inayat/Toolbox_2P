function for_behroo
file_name = 'all data_Behroo.xlsx';
[num,txt,raw] = xlsread(file_name,1,'A2:L26');

n = 0;
%%
num = [num num];
G1 = num(1:6,:);
G2 = num(7:13,:);
G3 = num(14:20,:);
G4 = num(21:25,:);

[within,dvn,xlabels,awithinD] = make_within_table({'S','B'},[4,6]);
dataT = make_between_table({G1;G2;G3;G4},dvn);
ra = RMA(dataT,within,{0.05,{'hsd','bonferroni'}});
ra.ranova

%%

G1 = num(1:6,1:6);
G2 = num(7:13,1:6);
G3 = num(14:20,1:6);
G4 = num(21:25,1:6);

[within,dvn,xlabels,awithinD] = make_within_table({'B'},[6]);
dataT = make_between_table({G1;G2;G3},dvn);
ra = RMA(dataT,within,{0.05,{'hsd','bonferroni'}});
ra.ranova
