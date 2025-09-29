function props = get_props_new(out,MV,props_i,si)

props = props_i;
n = 0;



% tCells = ((pvals < 0.05) * [4; 2; 1]) == 4; dCells = ((pvals < 0.05) * [4; 2; 1]) == 2; sCells = ((pvals < 0.05) * [4; 2; 1]) == 1;
% tdCells = ((pvals < 0.05) * [4; 2; 1]) == 6; dsCells = ((pvals < 0.05) * [4; 2; 1]) == 3; tsCells = ((pvals < 0.05) * [4; 2; 1]) == 5;
% ntdsCells = ((pvals < 0.05) * [4; 2; 1]) == 0; tdsCells = ((pvals < 0.05) * [4; 2; 1]) == 7;

ctnums = 0:3;
ctnums_names = {'NR','speed','time','time_speed'};
for cti = 1:length(ctnums)
    all_cells = {}; all_cells_zvals = {};
    for an = 1:5
        all_cells_an = {}; all_cells_an_zvals = {};
        for ii = 1:size(si,2)
            cn = si(1,ii); ap = si(2,ii);
            idx = 3; pvals = [MV{1}{an,cn,ap}.PC(:,idx) MV{2}{an,cn,ap}.PC(:,idx)];
            idx = 2; zvals = [MV{1}{an,cn,ap}.PC(:,idx) MV{2}{an,cn,ap}.PC(:,idx)];
            this_cells = ((pvals < 0.05) * [2; 1]) == ctnums(cti);
            all_cells_an = [all_cells_an this_cells];
            all_cells_an_zvals = [all_cells_an_zvals zvals];
        end
        all_cells = [all_cells;all_cells_an];
        all_cells_zvals = [all_cells_zvals;all_cells_an_zvals];
    end
    all_cells = cell_list_op(all_cells,props_i.good_FR,'and'); 
    cmdTxt = sprintf('props.newPC.cells_%s = all_cells;',ctnums_names{cti});
    eval(cmdTxt);
    cmdTxt = sprintf('props.newPC.cells_%s_zvals = all_cells_zvals;',ctnums_names{cti});
    eval(cmdTxt);

    all_cells = {}; all_cells_zvals = {};
    for an = 1:5
        all_cells_an = {}; all_cells_an_zvals = {};
        for ii = 1:size(si,2)
            cn = si(1,ii); ap = si(2,ii);
            idx = 3; pvals = [MV{1}{an,cn,ap}.MI(:,idx) MV{2}{an,cn,ap}.MI(:,idx)];
            idx = 2; zvals = [MV{1}{an,cn,ap}.MI(:,idx) MV{2}{an,cn,ap}.MI(:,idx)];
            this_cells = ((pvals < 0.05) * [2; 1]) == ctnums(cti);
            all_cells_an = [all_cells_an this_cells];
            all_cells_an_zvals = [all_cells_an_zvals zvals];
        end
        all_cells = [all_cells;all_cells_an];
        all_cells_zvals = [all_cells_zvals;all_cells_an_zvals];
    end
    all_cells = cell_list_op(all_cells,props_i.good_FR,'and'); 
    cmdTxt = sprintf('props.newMI.cells_%s = all_cells;',ctnums_names{cti});
    eval(cmdTxt);
    cmdTxt = sprintf('props.newMI.cells_%s_zvals = all_cells_zvals;',ctnums_names{cti});
    eval(cmdTxt);
end

cti = 2;
cmdTxt = sprintf('cells_or = props.newMI.cells_%s;',ctnums_names{cti}); eval(cmdTxt);
for cti = 3:4
    cmdTxt = sprintf('tempCe = props.newMI.cells_%s;',ctnums_names{cti}); eval(cmdTxt);
    cells_or = cell_list_op(cells_or,tempCe,'or');
end
props.newMI.cells_R = cells_or;

cti = 2;
cmdTxt = sprintf('cells_or = props.newPC.cells_%s;',ctnums_names{cti}); eval(cmdTxt);
for cti = 3:4
    cmdTxt = sprintf('tempCe = props.newPC.cells_%s;',ctnums_names{cti}); eval(cmdTxt);
    cells_or = cell_list_op(cells_or,tempCe,'or');
end
props.newPC.cells_R = cells_or;

% fields = fieldnames(props.newPC);

% for ii = 1:length(fields)
%     this_field = fields{ii};
%     cmdTxt = sprintf('tempCe = props.newPC.cells_%s;',this_field); eval(cmdTxt);
%     cmdTxt = sprintf('tempCe = props.newPC.cells_%s;',this_field); eval(cmdTxt);
% end
tProps = props;
singularly_tuned = cell_list_op(tProps.newPC.cells_time,tProps.newPC.cells_speed,'or');
mixed_tuned = tProps.newPC.cells_time_speed;
all_tuned = tProps.newPC.cells_time_speed;

props.newPC.cells_singularly_tuned = singularly_tuned;
props.newPC.cells_mixed_tuned = mixed_tuned;
props.newPC.cells_all_tuned = all_tuned;


singularly_tuned = cell_list_op(tProps.newMI.cells_time,tProps.newMI.cells_speed,'or');
mixed_tuned = tProps.newMI.cells_time_speed;
all_tuned = tProps.newMI.cells_time_speed;


props.newMI.cells_singularly_tuned = singularly_tuned;
props.newMI.cells_mixed_tuned = mixed_tuned;
props.newMI.cells_all_tuned = all_tuned;

props.newPC.names = {'NR','speed','time','time_speed','singularly_tuned','mixed_tuned','all_tuned'};
props.newMI.names = {'NR','speed','time','time_speed','singularly_tuned','mixed_tuned','all_tuned'};