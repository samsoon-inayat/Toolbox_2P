%% first find the variables time_cells, distance_cells, and speed_cells in the final_time_dist_bin...m file
variable_combs = {'FR_time','FR_dist','FR_speed'};
metric = 'MI';
metric = 'PC';
clear time_cells distance_cells speed_cells
clear t_cells_T d_cells_T s_cells_T td_cells_T ds_cells_T ts_cells_T tds_cells_T ntds_cells_T
clear t_cells_I d_cells_I s_cells_I td_cells_I ds_cells_I ts_cells_I tds_cells_I ntds_cells_I
all_cells = {};
ntrials = 30;
si = [Ar_t_T Ar_i_T ArL_t_T ArL_i_T Ars_t_T Ars_i_T]; 
Rs = o.Rs(:,si);
mRs = o.mR(:,si);

propsR = get_props_Rs(Rs,30);

mxsz = [];
for an = 1:5
    for ii = 1:length(si)
        tR = Rs{an,ii};
        mxsz(an,ii) = size(tR.sp_rasters,2);
    end
end

for an = 1:5
    for ii = 1:2:length(si)
        tR_on = Rs{an,ii};
        tR_off = Rs{an,ii+1};
        sp_rasters_on = tR_on.sp_rasters;
        sp_rasters_off = tR_off.sp_rasters;
        ncells = size(sp_rasters_off,3);
        for tr = 1:10
            nan_loc_on(tr) = find(isnan(squeeze(sp_rasters_on(tr,:,1))),1,'first');
            nan_loc_off(tr) = find(isnan(squeeze(sp_rasters_off(tr,:,1))),1,'first');
            if nan_loc_on(tr) > nan_loc_off(tr)
                the_lim(tr) = nan_loc_off(tr) - 1;
            else
                the_lim(tr) = nan_loc_on(tr) - 1;
            end
        end
        sp_rasters_on_lim = nan(10,max(the_lim),ncells);
        sp_rasters_off_lim = sp_rasters_on_lim;
        for tr = 1:10
            for cni = 1:ncells
                sp_rasters_on_lim(tr,1:the_lim(tr),cni) = sp_rasters_on(tr,1:the_lim(tr),cni);
                sp_rasters_off_lim(tr,1:the_lim(tr),cni) = sp_rasters_off(tr,1:the_lim(tr),cni);
            end
        end
        for cni = 1:ncells
            tR_on = squeeze(sp_rasters_on_lim(:,:,cni));
            tR_off = squeeze(sp_rasters_off_lim(:,:,cni));
            good_RF_on(cni) = sum(nansum(tR_on,2) > 0) > 3;
            good_RF_off(cni) = sum(nansum(tR_off,2) > 0) > 3;
        end
        perc_gRF(an,ii) = 100*sum(good_RF_on)/ncells;
        perc_gRF(an,ii+1) = 100*sum(good_RF_off)/ncells;
    end
end

avar = perc_gRF;
fac_names = {'CN','AP'}; fac_levels = [3,2];
[within,dvn,xlabels,awithinD] = make_within_table(fac_names,fac_levels);
dataT = make_between_table({avar},dvn);
ra = RMA(dataT,within,{0.05,{'hsd'}});
clc
print_for_manuscript(ra)


%%
si = [Ar_t_D Ar_i_D ArL_t_D ArL_i_D Ars_t_D Ars_i_D]; 
Rs = o.Rs(:,si);
mRs = o.mR(:,si);
for an = 1:5
    for ii = 1:length(si)
        tR = Rs{an,ii};
        for tii = 1:10
            
        end
    end
end