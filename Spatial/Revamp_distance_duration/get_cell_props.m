function get_cell_props(out,MV,var,selCells)

variable_combs = {'FR_time','FR_dist','FR_speed'};
vn = find(strcmp(var,variable_combs));
% Factors, CN, 
avar = [];
for an = 1:5
    anvar = [];
    for cn = 1:3
        for ap = 1:2
            AC = out.active_cells_trials{an,cn,ap};
            AAC1 = sum(AC,2) > 0 & sum(AC,2) <= 3;
            AAC2 = sum(AC,2) > 3 & sum(AC,2) <= 6;
            AAC3 = sum(AC,2) > 6 & sum(AC,2) < 10;
            selCells = AAC2 | AAC3; cellP = []; cellP1 = [];
            selCells = sum(AC,2) > 3;
            for vn = 1:length(variable_combs)
                tMV = MV{vn};
                tmet = tMV{an,cn,ap}.PC(selCells,3);
                tmet1 = tMV{an,cn,ap}.MI(selCells,3);
            %     tmet1 = tMV{an,cn,ap}.MI(AAC1,1); tmet2 = tMV{an,cn,ap}.MI(AAC2,1); tmet3 = tMV{an,cn,ap}.MI(AAC3,1);
            %     anvar = [anvar mean(tmet1) mean(tmet2) mean(tmet3)];
                cellP = [cellP tmet];
                cellP1 = [cellP1 tmet1];
            end
            cellT = cellP < 0.05;
            pTC = 100*(sum(cellT(:,1) & ~cellT(:,2) & ~cellT(:,3))/size(AC,1));
            pDC = 100*(sum(~cellT(:,1) & cellT(:,2) & ~cellT(:,3))/size(AC,1));
            pSC = 100*(sum(~cellT(:,1) & ~cellT(:,2) & cellT(:,3))/size(AC,1));
            pTDC = 100*(sum((cellT(:,1) | cellT(:,2)) & ~cellT(:,3))/size(AC,1));
            pTSC = 100*(sum((cellT(:,1) | cellT(:,3)) & ~cellT(:,2))/size(AC,1));
            pDSC = 100*(sum((cellT(:,2) | cellT(:,3)) & ~cellT(:,1))/size(AC,1));
            pTDS = 100*(sum((cellT(:,2) & cellT(:,3) & cellT(:,1)))/size(AC,1));
            anvar = [anvar pTC pDC pSC pTDC pTSC pDSC pTDS];

            cellT = cellP1 < 0.05;
            pTC = 100*(sum(cellT(:,1) & ~cellT(:,2) & ~cellT(:,3))/size(AC,1));
            pDC = 100*(sum(~cellT(:,1) & cellT(:,2) & ~cellT(:,3))/size(AC,1));
            pSC = 100*(sum(~cellT(:,1) & ~cellT(:,2) & cellT(:,3))/size(AC,1));
            pTDC = 100*(sum((cellT(:,1) | cellT(:,2)) & ~cellT(:,3))/size(AC,1));
            pTSC = 100*(sum((cellT(:,1) | cellT(:,3)) & ~cellT(:,2))/size(AC,1));
            pDSC = 100*(sum((cellT(:,2) | cellT(:,3)) & ~cellT(:,1))/size(AC,1));
            pTDS = 100*(sum((cellT(:,2) & cellT(:,3) & cellT(:,1)))/size(AC,1));
            anvar = [anvar pTC pDC pSC pTDC pTSC pDSC pTDS];
        end
    end
    avar = [avar;anvar];
end