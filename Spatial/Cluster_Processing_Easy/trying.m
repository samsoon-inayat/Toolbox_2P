function trying

%%
ii = 1;
onsets = ei{ii}.b.air_puff_r;
offsets = ei{ii}.b.air_puff_f;
dists = ei{ii}.b.dist(offsets) - ei{ii}.b.dist(onsets);