function showCellPopulations
%%
while 1
    an = 1; pl = 1;
    ntrials = 50; %si = [Lb_T ArL_L_T Lbs_T Ab_t_T Ab_i_T Abs_t_T Abs_i_T Ar_t_D ArL_t_D Ars_t_D Ar_t_T ArL_t_T Ars_t_T Ar_i_D ArL_i_D Ars_i_D Ar_i_T ArL_i_T Ars_i_T];
    si = [Lb_T ArL_L_T Lbs_T Ab_T Abs_T Ar_t_D ArL_t_D Ars_t_D Ar_i_T ArL_i_T Ars_i_T];
    tRs = o.Rs(:,si);
    props1 = get_props_Rs(tRs,ntrials);
    resp = props1.good_FR;

    tei = ei{an};
    reD = get_cell_list(resp,[6 7 8]);
    reT = get_cell_list(resp,[9 10 11]);
    for ii = 1:size(resp,2)
        selCellsD{ii} = get_cell_nums(tRs{an,ii},reD{an,ii},pl);
        selCellsT{ii} = get_cell_nums(tRs{an,ii},reT{an,ii},pl);
        selCells{ii} = get_cell_nums(tRs{an,ii},resp{an,ii},pl);
    end
    
    break;
end

figure(1000);clf;
ha = axes;
showCells(ha,tei,pl,selCells([6 7 8]),[0.3 0.3]);
% showCells(ha,tei,pl,{selCellsD{1} selCellsT{1}},[0.3 0.3]);


