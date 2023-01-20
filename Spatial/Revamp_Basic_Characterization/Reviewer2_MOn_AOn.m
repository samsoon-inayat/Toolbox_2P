%% Three-Way RM-ANOVA Most Recent compare all variables pooled air on or off brake vs no-brake exc and inh (and controls) big ANOVA test

    ntrials = 50; %si = [Lb_T ArL_L_T Lbs_T Ab_t_T Ab_i_T Abs_t_T Abs_i_T Ar_t_D ArL_t_D Ars_t_D Ar_t_T ArL_t_T Ars_t_T Ar_i_D ArL_i_D Ars_i_D Ar_i_T ArL_i_T Ars_i_T];
    event_type = {'2-B-AOn','7-B-AOn','2-B-AOff','7-B-AOff','2-B-Arb','7-B-Arb','3-NB-AOn','4-NB-AOn','5-NB-AOn','3-NB-AOff',...
        '4-NB-AOff','5-NB-AOff','3-NB-Arb','4-NB-Arb','5-NB-Arb','1-B-L','6-B-L','3-NB-A','5-NB-A','4-NB-AL'};
    event_type = {'B-AOn-Exc','B-AOn-Inh','B-AOff-Exc','B-AOff-Inh','B-Arb-Exc','B-Arb-Inh','NB-AOn-Exc','NB-AOn-Inh','NB-AOff-Exc','NB-AOff-Inh','NB-Arb-Exc','NB-Arb-Inh'};
    sic = {[Ab_On Abs_On];[Ab_Off Abs_Off];[Ab_Offc Abs_Offc];[Ar_On ArL_On Ars_On];[Ar_Off ArL_Off Ars_Off];[Ar_Offc ArL_Offc Ars_Offc]};
    sic = {[Ar_On ArL_On Ars_On];[M_Off]};
%     sic = {Ar_On;ArL_On;Ars_On;M_On};
%     sic = {[Ab_On];[Abs_On];[Ab_Off];[Abs_Off];[Ab_Offc];[Abs_Offc];[Ar_On];[ArL_On];[Ars_On];[Ar_Off];[ArL_Off];[Ars_Off];[Ar_Offc];[ArL_Offc];[Ars_Offc];[Lb];[Lbs];[Ar_L];[Ars_L];[Ar_L]}; % for heat map after light sitmulus
    
%     event_type = {'2-B-AOn','7-B-AOn','2-B-AOff','7-B-AOff','3-NB-AOn','4-NB-AOn','5-NB-AOn','3-NB-AOff',...
%         '4-NB-AOff','5-NB-AOff','1-B-L','6-B-L','4-NB-AL','3-NB-A','5-NB-A'};
%     sic = {[Ab_On];[Abs_On];[Ab_Off];[Abs_Off];[Ar_On];[ArL_On];[Ars_On];[Ar_Off];[ArL_Off];[Ars_Off];[Lb];[Lbs];[ArL_L];[Ar_L];[Ars_L]}; % for heat map after light sitmulus
%     sic = {[Ab_On Abs_On];[Ar_On ArL_On Ars_On];[Ab_Off Abs_Off];[Ar_Off ArL_Off Ars_Off];[Ab_Offc Abs_Offc];[Ar_Offc ArL_Offc Ars_Offc]};
    clear all_gFR all_exc all_inh all_gV all_exc_inh
    prop_names = {'resp','N_Resp_Trials','zMI','zMINaN','HaFD','HiFD','cells_pooled'};
    cell_sel = {'good_FR','vals','exc','inh'};
    varName = {'all_gFR','all_gV','all_exc','all_inh'};
    for cii = 1:length(cell_sel);
        cmdTxt = sprintf('clear %s',varName{cii});eval(cmdTxt);
    end
    pni = 7;
    all_exc_inh = [];
    for ii = 1:length(sic)
        sit = sic{ii};
        tRs = o.Rs(:,sit);
        props1 = get_props_Rs(tRs,ntrials);
        for cii = 1:length(cell_sel)
            cmdTxt = sprintf('gFR = props1.%s;',cell_sel{cii});eval(cmdTxt);
            if pni == 1
                rf = find_percent(gFR);
                cmdTxt = sprintf('%s(:,ii) = mean(rf,2);',varName{cii}); eval(cmdTxt)
            else
                if pni == 7
                    cmdTxt = sprintf('%s(:,ii) = cell_list_op(gFR,[],''or'',1);',varName{cii});eval(cmdTxt);
                else
                    if pni == 3 || pni == 4
                        cmdtxt = sprintf('rf = props1.%s;',prop_names{3}); eval(cmdtxt);
                        if pni == 3
                            cmdTxt = sprintf('%s(:,ii) = mean(exec_fun_on_cell_mat(rf,''nanmean'',gFR),2);',varName{cii});eval(cmdTxt);
                        else
                            for rrr = 1:size(rf,1)
                                for ccc = 1:size(rf,2)
                                    temp = isnan(rf{rrr,ccc});
                                    rf1(rrr,ccc) = 100*sum(temp(gFR{rrr,ccc}))/size(rf{rrr,ccc},1);
                                end
                            end
                            cmdTxt = sprintf('%s(:,ii) = mean(rf1,2);',varName{cii}); eval(cmdTxt)
                        end
                    else
                        cmdtxt = sprintf('rf = props1.%s;',prop_names{pni}); eval(cmdtxt);
                        cmdTxt = sprintf('%s(:,ii) = mean(exec_fun_on_cell_mat(rf,''mean'',gFR),2);',varName{cii});eval(cmdTxt);
                    end
                end
            end
        end
        all_exc_inh = [all_exc_inh all_exc(:,ii) all_inh(:,ii)];
    end
    %%
     %
    avar = find_percent(all_gV);
    [within,dvn,xlabels,withinD] = make_within_table({'Cond'},[2]); withinD3 = withinD;
    dataT = make_between_table({avar},dvn);
    ra = RMA(dataT,within,{0.05,{'hsd'}});
    ra.ranova
    print_for_manuscript(ra)
    
    [h,p,cd,stat] = ttest(avar(:,1),avar(:,2))
    %%
    