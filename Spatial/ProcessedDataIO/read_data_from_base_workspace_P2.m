function out = read_data_from_base_workspace_P2(selC)

rsel = 1
ei_names = {'ei'};
ET_names = {'ET'};
param_names = {'P2_1'};

for ii = 1:length(ei_names)
    rsel = ii;
    all_ETs{ii} = evalin('base',ET_names{ii});
    all_eis{ii} = evalin('base',sprintf('%s',ei_names{rsel}));
    all_paramMs{ii} = parameter_matrices_ctrl('get',sprintf('%s',param_names{rsel}));
%     all_paramMs{ii}.belt_lengths = get_belt_length(all_eis{ii})
    [all_cpMs{ii},all_pMs{ii}] = parameter_matrices_ctrl('select',sprintf('%s',param_names{rsel}),{all_paramMs{ii},selC});
    all_selAnimals{ii} = 1:length(all_eis{ii})
    all_perc_cells{ii} = parameter_matrices_ctrl('print',sprintf('%s',param_names{rsel}),{all_cpMs{ii},all_pMs{ii},all_ETs{ii},all_selAnimals{ii}});
end

out.ei_names = ei_names; out.ET_names = ET_names; out.param_names = param_names;
out.ETs = all_ETs; out.eis = all_eis; out.paramMs = all_paramMs;
out.cpMs = all_cpMs; out.pMs = all_pMs;
out.selAnimals = all_selAnimals;
out.perc_cells = all_perc_cells;

