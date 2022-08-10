% disregard this code for now Aug 2022

%% general for all properties including responsivity, response fidelity, zMI, Rs

    ntrials = 50;
    si = {[Ab_On Abs_On];[Ab_Off Abs_Off];[Ab_Offc Abs_Offc];[Ar_On ArL_On Ars_On];[Ar_Off ArL_Off Ars_Off];[Ar_Offc ArL_Offc Ars_Offc]};
%     si = [C1_i_T C2_i_T C3_i_T C4_i_T];
    Rs_C = oC.Rs(:,si); Rs_A = oA.Rs(:,si); mRs_C = oC.mR(:,si); mRs_A = oA.mR(:,si);
    props_C = get_props_Rs(oC.Rs(:,si),ntrials); props_A = get_props_Rs(oA.Rs(:,si),ntrials);
    pop_var_name = {'all','vals','valsT','Nvals','good_zMI','Ngood_zMI'};
    pop_var_name = {'good_zMI','good_Gauss','good_MFR'};
    pop_var_name = {'all'};
    pop_var_name = {'vals'};
    sel_pop_C = cell_list_op(props_C,pop_var_name); sel_pop_A = cell_list_op(props_A,pop_var_name);
    cell_types = {'C1','C2','C3','C4'};
    
    params = {'perc','N_Resp_Trials','zMI','rs','nan_zMI','nan_rs','HaFD','HiFD','PWs','centers','peak_locations','mean_FR','MFR'};
    varT = 3;%:length(params)
    for pii = varT
        if pii == 1
            mean_var_C = exec_fun_on_cell_mat(sel_pop_C,'percent'); mean_var_A = exec_fun_on_cell_mat(sel_pop_A,'percent'); 
        else
            eval(sprintf('var_C = get_vals(props_C.%s,sel_pop_C);',params{pii})); eval(sprintf('var_A = get_vals(props_A.%s,sel_pop_A);',params{pii}));
            if pii == 5 || pii == 6
                mean_var_C = exec_fun_on_cell_mat(sel_pop_C,'percent'); mean_var_A = exec_fun_on_cell_mat(sel_pop_A,'percent'); 
            else
                mean_var_C = exec_fun_on_cell_mat(var_C,'nanmean'); mean_var_A = exec_fun_on_cell_mat(var_A,'nanmean'); 
            end
        end
    end
    varC = mean_var_C;
    varA = mean_var_A;
    [within,dvn,xlabels] = make_within_table({'Cond'},[4]);
    dataT = make_between_table({varC;varA},dvn);
    ra = RMA(dataT,within,{0.05,{'bonferroni'}});
    ra.ranova
