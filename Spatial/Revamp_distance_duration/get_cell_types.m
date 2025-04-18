function ctout = get_cell_types(out,MV)

variable_combs = {'FR_time','FR_dist','FR_speed'};
cell_types = {'NoT','speed','dist','dist-speed','time','time-speed','time-dist','TDS'};
ctout.cell_types = cell_types;
for an = 1:5
    for cn = 1:3
        for ap = 1:2
            AC = out.active_cells_trials{an,cn,ap};
            ctout.RF1{an,cn,ap} = sum(AC,2) > 1 & sum(AC,2) <= 4;
            ctout.RF2{an,cn,ap}  = sum(AC,2) > 4 & sum(AC,2) <= 7;
            ctout.RF3{an,cn,ap}  = sum(AC,2) > 7 & sum(AC,2) < 10;
            ctout.RFg5{an,cn,ap}  = sum(AC,2) > 5;
            ctout.RFl5{an,cn,ap}  = sum(AC,2) <= 5;
            cellP = []; cellP1 =[];
            for vn = 1:length(variable_combs)
                tMV = MV{vn};
                tmet = tMV{an,cn,ap}.PC(:,3);
                tmet1 = tMV{an,cn,ap}.MI(:,3);
                cellP = [cellP tmet];
                cellP1 = [cellP1 tmet1];
            end
            cellT = cellP < 0.05;
            dcpc = cellT * [4; 2; 1];
            ctout.CT_PC{an,cn,ap} = cell_types(dcpc+1);
            ctout.CT_PCd{an,cn,ap} = dcpc+1;

            cellT = cellP1 < 0.05;
            dcmi = cellT * [4; 2; 1];
            ctout.CT_MI{an,cn,ap} = cell_types(dcmi+1);
            ctout.CT_MId{an,cn,ap} = dcmi+1;
        end
    end
end
