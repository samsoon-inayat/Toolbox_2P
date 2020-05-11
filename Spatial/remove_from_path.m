function remove_from_path
selpp = []; selpf = [];

pp_folder_names = {'PostProcessing','PostProcessing_15','PostProcessing_16'};
pf_folder_names = {'Protocol10_Figures','Protocol10_Figures_AD','Protocol15_1_Figures','Protocol15_Figures','Protocol_16_Figures'};

rmp(pp_folder_names,selpp);
rmp(pf_folder_names,selpf)

function rmp(pp_folder_names,selpp)
tot = 1:length(pp_folder_names);
remaining = setdiff(tot,selpp);
for ii = 1:length(remaining)
    p = genpath(fullfile(pwd,pp_folder_names{remaining(ii)}));
    rmpath(p);
end


