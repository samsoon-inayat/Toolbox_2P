function figure_place_cells_raster_types_conjuctions(varargin)

owr = 0;
ei = evalin('base','ei10');
selAnimals = [1:4 9];
% in the following variable all the measurements are in the matrices form
% for each variable colums indicate raster and stim marker types specified 
% the rows indicate condition numbers.
paramMs = get_parameters_matrices(ei,selAnimals,owr);
% after getting all matrics, we can apply selection criteria to select a
% subgroup of cells
% here is the selection criteria

selC.areCells = 1; % to see if it is identified as cell by Suite2P
selC.plane_number = NaN; % to select a particular plane to analyze
selC.conditionNumber = [2];
selC.rasterType = [1];
selC.zMI_threshold = 3; % to select tuned cells
selC.fwidth_limits = NaN; % to select limits of field widths to be included
selC.fcenter_limits = NaN; % to select field center locations to be included
selC.frs_threshold = NaN; % to select a threshold to be put on fitting r-square value


paramMs1 = get_parameters_matrices(paramMs,selC);







mData = evalin('base','mData');
colors = mData.colors;
sigColor = mData.sigColor;

mData.belt_length = ei{selAnimals(1)}.b.belt_length;
n = 0;

%%
% there are four raster types
% there are four condition types
% for each raster type there are four conditions
% in total therefore are 16 combinations
% from the point of view of a cell, there are four raster types and four
% conditions. Therefore I have a 4x4 matrix showing zMI and pc vals based
% on a threshold of zMI = 3
% thus for each animal I will have a 4x4xNcells matrix

for an = 1:length(selAnimals)
    pc_M = []; zMI_M = [];
    for ii = 1:size(all_pcs,2) % correspond to condition numbers
        all_pcs_c = all_pcs{ii};
        all_zMIs_c = all_zMIs{ii};
        one_an = all_pcs_c(an,:);
        one_an_z = all_zMIs_c(an,:);
        for jj = 1:size(one_an,2) % correspond to 4 raster types
            pc_M(ii,jj,:) = one_an{jj};
            zMI_M(ii,jj,:) = one_an_z{jj};
        end
    end
    pc{an} = pc_M;
    zMI{an} = zMI_M;
end

%%
mdthiszMI = [];
for an = 1:length(pc)
    thispc = pc{an};
    thiszMI = zMI{an};
    dthiszMI_c = diff(thiszMI,1);
    mdthiszMI(:,:,an) = nanmean(dthiszMI_c,3);
    rthiszMI = reshape(thiszMI,size(thiszMI,1)*size(thiszMI,2),size(thiszMI,3))';
end

[avg,sem] = findMeanAndStandardError(mdthiszMI,3);