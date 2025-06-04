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
%% Find time cells with PC in air_on phase and plot graphs and look at population stats
% this is with time-binning
% RV is response variable
tt = 1;
clear pvals zvals mi_zvals tCells
RV = 'perc';
% RV = 'rfid';
% RV = 'mval';
% RV = 'MI';
% RV = 'PC';
avar = [];
for an = 1:5
    anvar = [];
    for cn = 1:3
        ap = 1;
        AC = outT.active_cells_trials{an,cn,ap}; mAC1 = outT.active_cells{an,cn,ap}; RF = sum(AC,2); mAC = RF>5;
        out = outT; MV = met_valsT;         % out = outD; MV = met_valsD; 
        selcells = mAC1;% & (RF>0);
        idx = 3; pvals = [MV{tt}{an,cn,ap}.PC(:,idx)]; idx = 2; zvals = [MV{tt}{an,cn,ap}.PC(:,idx)];
        idx = 2; mi_zvals = [MV{tt}{an,cn,ap}.MI(:,idx)];
        % idx = 3; pvals = [MV{1}{an,cn,ap}.MI(:,idx)];
        % idx = 2; zvals = [MV{1}{an,cn,ap}.MI(:,idx)];
        % criteria for finding the types of cells
        tCells = pvals < 0.05; %time_cells{an,cn,ap,bn} = tCells< 0.05;
        if strcmp(RV,'perc')
            anvar = [anvar sum(tCells & selcells)/length(AC)];
        end
        if strcmp(RV,'mval')
            anvar = [anvar nanmean(zvals(tCells&selcells,1))];
        end
        if strcmp(RV,'rfid')
            anvar = [anvar mean(RF(tCells&selcells,1))];
        end
        if strcmp(RV,'MI')
            anvar = [anvar mean(mi_zvals(tCells&selcells,1))];
        end
        if sum(isnan(anvar)) > 0
            [an cn]
        end
    end
    if strcmp(RV,'perc')
        avar = [avar;100*anvar];
    else
        avar = [avar;anvar];
    end
end

fac_names = {'CN'}; fac_levels = [3];
% fac_names = {'CN','AP','BT','RF','CT'}; fac_levels = [3,2,2,3,4];
% fac_names = {'CN','AP','BT','CT'}; fac_levels = [3,2,2,7];
[within,dvn,xlabels,awithinD] = make_within_table(fac_names,fac_levels);
dataT = make_between_table({avar},dvn);
ra = RMA(dataT,within,{0.05,{''}});
clc
print_for_manuscript(ra)
% (Intercept) [F(1,4) = 93.20, p < 0.001, η2 = .90] <--
%% visualizing the results in the previous section
tcolors = repmat(mData.dcolors(8:10),1,3); MY = 30; ysp = 1; mY = 0; ystf = 2; ysigf = 0.025;titletxt = ''; ylabeltxt = {'Cells (%)'}; % for all cells (vals) MY = 80
[hbs,xdata,mVar,semVar,combs,p,h] = view_results_rmanova([],ra,{'CN','hsd',0.05},[1 1.75],tcolors,[mY MY ysp ystf ysigf],mData);