%% first find the variables time_cells, distance_cells, and speed_cells in the final_time_dist_bin...m file
variable_combs = {'FR_time','FR_dist','FR_speed'};
% metric = 'MI';
metric = 'PC';
for an = 1:5
    for cn = 1:3
        for ap = 1:2
            AC = outT.active_cells_trials{an,cn,ap}; mAC = outT.active_cells{an,cn,ap}; RF = sum(AC,2);
            for bn = 1:2
                if bn == 1
                    out = outT; MV = met_valsT; 
                else
                    out = outD; MV = met_valsD; 
                end
                selcells = mAC;% & (RF>0);
                if strcmp(metric,'PC')
                    idx = 3; pvals = [MV{1}{an,cn,ap}.PC(:,idx) MV{2}{an,cn,ap}.PC(:,idx) MV{3}{an,cn,ap}.PC(:,idx)];
                    idx = 2; zvals = [MV{1}{an,cn,ap}.PC(:,idx) MV{2}{an,cn,ap}.PC(:,idx) MV{3}{an,cn,ap}.PC(:,idx)];
                else
                    idx = 3; pvals = [MV{1}{an,cn,ap}.MI(:,idx) MV{2}{an,cn,ap}.MI(:,idx) MV{3}{an,cn,ap}.MI(:,idx)];
                    idx = 2; zvals = [MV{1}{an,cn,ap}.MI(:,idx) MV{2}{an,cn,ap}.MI(:,idx) MV{3}{an,cn,ap}.MI(:,idx)];
                end
                % criteria for finding the types of cells
                tCells = ((pvals < 0.05) * [4; 2; 1]) == 4; dCells = ((pvals < 0.05) * [4; 2; 1]) == 2; sCells = ((pvals < 0.05) * [4; 2; 1]) == 1;
                time_cells{an,cn,ap,bn} = tCells; distance_cells{an,cn,ap,bn} = dCells; speed_cells{an,cn,ap,bn} = sCells;
            end
        end
    end
end
%% Plot rasters for one animal top-down concatenated time-distance-and speed cells
% out = outT; bn = 1;
out = outD; bn = 2;
for an = 1%:5
    for cn = 3%:3
        for ap = 1
            rasters = out.arasters{an,cn,ap}; rx = out.arx{an,cn,ap}; ry = out.ary{an,cn,ap};
            tCells = time_cells{an,cn,ap,bn};
            dCells = distance_cells{an,cn,ap,bn};
            sCells = speed_cells{an,cn,ap,bn};
            num_neurons = size(rasters,3);
            for nn = 1:num_neurons
                if tCells(nn) == 0
                    continue;
                end
                figure(100);clf;
                imagesc(rx,ry,rasters(:,:,nn));colorbar;
                title(nn)
                pause(0.3);
            end
        end
    end
end

%%