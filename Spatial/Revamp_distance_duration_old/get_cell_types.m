function ctout = get_cell_types(out,MV,alpha)

variable_combs = {'FR_time','FR_dist','FR_speed'};
cell_types = {'NoT','speed','dist','dist-speed','time','time-speed','time-dist','TDS'};
ctout.cell_types = cell_types;
for an = 1:5
    for cn = 1:3
        for ap = 1:2
            AC = out.active_cells_trials{an,cn,ap};
            ctout.RF1{an,cn,ap} = sum(AC,2) > 0 & sum(AC,2) <= 4;
            ctout.RF2{an,cn,ap}  = sum(AC,2) > 4 & sum(AC,2) <= 7;
            ctout.RF3{an,cn,ap}  = sum(AC,2) > 7 & sum(AC,2) <= 10;
            cellP = []; cellP1 =[]; cellP2 =[]; cellP3 =[]; all_cellP =[];
            for vn = 1:length(variable_combs)
                tMV = MV{vn};
                tmet = tMV{an,cn,ap}.PC(:,3);
                tmet1 = tMV{an,cn,ap}.MI(:,3);
                cellP = [cellP tmet];
                cellP1 = [cellP1 tmet1];
                tmet2 = tMV{an,cn,ap}.PC(:,2);
                tmet3 = tMV{an,cn,ap}.MI(:,2);
                cellP2 = [cellP2 tmet2];
                cellP3 = [cellP3 tmet3];
            end
            cellT = cellP < alpha;
            % ctout.CT_Tt{an,cn,ap} = get_tuning_type(cellP,cellP2,cellP3);
            dcpc = cellT * [4; 2; 1];
            ctout.CT_PC{an,cn,ap} = cell_types(dcpc+1);
            ctout.CT_PCd{an,cn,ap} = dcpc+1;

            cellT = cellP1 < alpha;
            dcmi = cellT * [4; 2; 1];
            ctout.CT_MI{an,cn,ap} = cell_types(dcmi+1);
            ctout.CT_MId{an,cn,ap} = dcmi+1;

            cellT = cellP1 < alpha & cellP < alpha;
            dcboth = cellT * [4; 2; 1];
            ctout.CT_B{an,cn,ap} = cell_types(dcboth+1);
            ctout.CT_Bd{an,cn,ap} = dcmi+1;
        end
    end
end


function Tt = get_tuning_type(pvals,pcs,mis)

n = 0;

plogical = pvals < 0.05;
dcpc = plogical * [4; 2; 1];
mixtuning = dcpc > 1;

err



% (Intercept) [F(1,4) = 151.20, p < 0.001, η2 = .69] <--
% (Intercept):BT [F(1,4) = 10.11, p = .034, η2 = .03] <--
% (Intercept):RF [F(1,4) = 101.48, p < 0.001, η2 = .33] <--
% (Intercept):CT [F(2,8) = 10.01, p = .007, η2 = .12] <--
% (Intercept):CN:BT [F(2,8) = 5.90, p = .027, η2 = .01] <--
% (Intercept):AP:RF [F(1,4) = 35.19, p = .004, η2 = .08] <--
% (Intercept):AP:CT [F(2,8) = 4.71, p = .045, η2 = .04] <--
% (Intercept):BT:CT [F(2,8) = 7.65, p = .014, η2 = .05] <--
% (Intercept):MT:CT [F(2,8) = 20.03, p < 0.001, η2 = .16] <--
% (Intercept):CN:BT:MT [F(2,8) = 8.66, p = .010, η2 = .01] <--
% (Intercept):AP:BT:MT [F(1,4) = 25.01, p = .007, η2 = .01] <--
% (Intercept):AP:BT:CT [F(2,8) = 9.61, p = .007, η2 = .02] <--
% (Intercept):BT:RF:CT [F(2,8) = 5.95, p = .026, η2 = .01] <--
% (Intercept):BT:MT:CT [F(2,8) = 8.24, p = .011, η2 = .03] <--
% (Intercept):AP:BT:RF:CT [F(2,8) = 8.51, p = .010, η2 = .02] <--
% (Intercept):AP:BT:MT:CT [F(2,8) = 6.10, p = .025, η2 = .01] <--
% (Intercept):BT:RF:MT:CT [F(2,8) = 11.33, p = .005, η2 = .00] <--