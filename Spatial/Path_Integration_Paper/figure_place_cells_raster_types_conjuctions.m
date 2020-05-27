function figure_place_cells_raster_types_conjuctions(varargin)

[all_pcs,all_zMIs] = find_all_pcs_zMIs;
ei = evalin('base','ei10');
mData = evalin('base','mData');
colors = mData.colors;
sigColor = mData.sigColor;
selAnimals = [1:9];
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