function processing_trial_by_trial
%%
% G = 'C'; o = oC;
% G = 'A'; o = oA;

si = [C1_t_D C1_i_T C2_t_D C2_i_T C3_t_D C3_i_T C4_t_D C4_i_T];
[allRsC_C,allmRsT_C,allresp_C,pcs_C] = return_values_trial_to_trial_Analysis(oC,si);

[allRsC_A,allmRsT_A,allresp_A,pcs_A] = return_values_trial_to_trial_Analysis(oA,si);

%% for responsivity
varC = pcs_C.respRV;
varA = pcs_A.respRV;
[within,dvn,xlabels] = make_within_table({'Cond','DT','Trials'},[4,2,10]);
dataT = make_between_table({varC;varA},dvn);
ra = RMA(dataT,within,{0.05,{'hsd','bonferroni'}});
ra.ranova


%% for conjunction

varC = []; varA = [];
for an = 1:5
    tvar = pcs_C.conjV(an,:);
    tvar(10:10:79) = NaN; tvar(isnan(tvar)) = []; 
    varC(an,:) = tvar;
    tvar = pcs_A.conjV(an,:);
    tvar(10:10:79) = NaN; tvar(isnan(tvar)) = []; 
    varA(an,:) = tvar;
end

[within,dvn,xlabels] = make_within_table({'Cond','DT','Trials'},[4,2,9]);
dataT = make_between_table({varC;varA},dvn);
ra = RMA(dataT,within,{0.05,{'hsd','bonferroni'}});
ra.ranova