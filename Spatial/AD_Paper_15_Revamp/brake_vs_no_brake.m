function brake_vs_no_brake
%%
ntrials = 50;
event_type = {'B-AOn-Exc','B-AOn-Inh','B-AOff-Exc','B-AOff-Inh','B-Arb-Exc','B-Arb-Inh','NB-AOn-Exc','NB-AOn-Inh','NB-AOff-Exc','NB-AOff-Inh','NB-Arb-Exc','NB-Arb-Inh'};
sic = {[Ab_On Abs_On];[Ab_Off Abs_Off];[Ab_Offc Abs_Offc];[Ar_On ArL_On Ars_On];[Ar_Off ArL_Off Ars_Off];[Ar_Offc ArL_Offc Ars_Offc]};


event_type = {'1-L','2-AOn','2-AOff','2-Arb','3-AOn','3-A','3-AOff','3-Arb','4-AOn','4-AL','4-AOff','4-Arb','5-AOn','5-A','5-AOff','5-Arb','6-L','7-AOn','7-AOff','7-Arb'};
sic = {Lb;Ab_On;Ab_Off;Ab_Offc;Ar_On;Ar_L;Ar_Off;Ar_Offc;ArL_On;ArL_L;ArL_Off;ArL_Offc;Ars_On;Ars_L;Ars_Off;Ars_Offc;Lbs;Abs_On;Abs_Off;Abs_Offc};

event_type = {'1-L','2-AOn','2-AOff','3-AOn','3-AOff','4-AOn','4-AL','4-AOff','5-AOn','5-AOff','6-L','7-AOn','7-AOff'};
sic = {Lb;Ab_On;Ab_Off;Ar_On;Ar_Off;ArL_On;ArL_L;ArL_Off;Ars_On;Ars_Off;Lbs;Abs_On;Abs_Off};

pni = 7;
[all_gFR_C,all_gV_C,all_exc_C,all_inh_C,all_exc_inh_C] = return_values_props(oC,sic,pni);
[all_gFR_A,all_gV_A,all_exc_A,all_inh_A,all_exc_inh_A] = return_values_props(oA,sic,pni);

%%
varC = all_exc_inh_C;
varA = all_exc_inh_A;
[within,dvn,xlabels,awithinD] = make_within_table({'ST','ET','CT'},[2,3,2]);
dataT = make_between_table({varC;varA},dvn);
ra = RMA(dataT,within,{0.05,{'hsd','bonferroni'}});
ra.ranova

%%
while 1
    varC = all_gV_C;
    varA = all_gV_A;
    [within,dvn,xlabels,awithinD] = make_within_table({'ST','ET'},[2,3]);
    dataT = make_between_table({varC;varA},dvn);
    ra = RMA(dataT,within,{0.05,{'hsd','bonferroni'}});
    ra.ranova
    %% 2nd level ET = 1
    varC1 = varC(:,awithinD(:,2)==1);
    varA1 = varA(:,awithinD(:,2)==1);
    [within,dvn,xlabels,withinD] = make_within_table({'ST'},[2]);
    dataT = make_between_table({varC1;varA1},dvn);
    ra1 = RMA(dataT,within,{0.05,{'hsd','bonferroni'}});
    ra1.ranova
    
    %% 2nd level ET = 2
    varC1 = varC(:,awithinD(:,2)==2);
    varA1 = varA(:,awithinD(:,2)==2);
    [within,dvn,xlabels,withinD] = make_within_table({'ST'},[2]);
    dataT = make_between_table({varC1;varA1},dvn);
    ra1 = RMA(dataT,within,{0.05,{'hsd','bonferroni'}});
    ra1.ranova
    
    %% 2nd level ET = 3
    varC1 = varC(:,awithinD(:,2)==3);
    varA1 = varA(:,awithinD(:,2)==3);
    [within,dvn,xlabels,withinD] = make_within_table({'ST'},[2]);
    dataT = make_between_table({varC1;varA1},dvn);
    ra1 = RMA(dataT,within,{0.05,{'hsd','bonferroni'}});
    ra1.ranova
    
    %% 2nd level ST = 1
    varC1 = varC(:,awithinD(:,1)==1);
    varA1 = varA(:,awithinD(:,1)==1);
    [within,dvn,xlabels,withinD] = make_within_table({'ET'},[3]);
    dataT = make_between_table({varC1;varA1},dvn);
    ra1 = RMA(dataT,within,{0.05,{'hsd','bonferroni'}});
    ra1.ranova
    %% 2nd level ST = 2
    varC1 = varC(:,awithinD(:,1)==2);
    varA1 = varA(:,awithinD(:,1)==2);
    [within,dvn,xlabels,withinD] = make_within_table({'ET'},[3]);
    dataT = make_between_table({varC1;varA1},dvn);
    ra1 = RMA(dataT,within,{0.05,{'hsd','bonferroni'}});
    ra1.ranova
    %% 2nd level Group = 1
    [within,dvn,xlabels,withinD] = make_within_table({'ST','ET'},[2,3]);
    dataT = make_between_table({varC},dvn);
    ra11 = RMA(dataT,within,{0.05,{'hsd','bonferroni'}});
    ra11.ranova
    
    [within,dvn,xlabels,withinD] = make_within_table({'ST','ET'},[2,3]);
    dataT = make_between_table({varA},dvn);
    ra12 = RMA(dataT,within,{0.05,{'hsd','bonferroni'}});
    ra12.ranova
%%
break
end
