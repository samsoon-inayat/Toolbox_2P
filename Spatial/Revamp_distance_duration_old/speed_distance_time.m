an = 1;

data_an = udataT{an};
redges = find_rising_edge(data_an.bnb,0.5,0.05);
fedges = find_falling_edge(data_an.bnb,-0.5,0.05);
data_an.nb = zeros(size(data_an.bnb));
data_an.nb((fedges(1)-1):(redges(2)+1))=1;
figure(100);clf;plot(data_an.ts,[data_an.nb',data_an.bnb']);pause(0.1)
%
r_data_an = data_an;
fields = fieldnames(data_an);
idx = data_an.nb == 1; % logical index
for ii = 1:numel(fields)
    if length(data_an.(fields{ii})) == length(idx)
        r_data_an.(fields{ii}) = data_an.(fields{ii})(:,idx);
    end
end
num_neurons = size(r_data_an.firing_rate,1);
ts = r_data_an.ts;ds = r_data_an.ds;sp = r_data_an.speed;tm = r_data_an.animal_motion;
figure(100);clf;plot(ts,r_data_an.air & r_data_an.C3)
% figure(100);clf;plot(ts,~r_data_an.air);hold on;
% plot(ts,r_data_an.air);
% %%
% C3_air_on = r_data_an.air & r_data_an.C3;
% figure(100);clf;plot(ts,C3_air_on)
% tts = ts(C3_air_on == 1); tds = ds(C3_air_on == 1); tsp = sp(C3_air_on == 1); ttm = tm(C3_air_on == 1);
% FR = r_data_an.firing_rate(:,C3_air_on == 1);
% out = find_cell_speed_tuning(FR,tsp);

%%
csel = r_data_an.C3;% & r_data_an.air;
figure(100);clf;plot(ts,csel);pause(0.5);
tts = ts(csel == 1); tds = ds(csel == 1); tsp = sp(csel == 1); ttm = tm(csel == 1);
FR = r_data_an.firing_rate(:,csel == 1);
out3 = find_cell_speed_tuning(FR,tsp);

csel = r_data_an.C4;% & r_data_an.air;
figure(100);clf;plot(ts,csel);pause(0.5);
tts = ts(csel == 1); tds = ds(csel == 1); tsp = sp(csel == 1); ttm = tm(csel == 1);
FR = r_data_an.firing_rate(:,csel == 1);
out4 = find_cell_speed_tuning(FR,tsp);

csel = r_data_an.C5;% & r_data_an.air;
figure(100);clf;plot(ts,csel);pause(0.5);
tts = ts(csel == 1); tds = ds(csel == 1); tsp = sp(csel == 1); ttm = tm(csel == 1);
FR = r_data_an.firing_rate(:,csel == 1);
out5 = find_cell_speed_tuning(FR,tsp);

%%
d.bcs = out3.bin_centers;
fitg = out3.fits.gauss; [rs,MFR,centers,PWs] = get_gauss_fit_parameters(fitg.coeffsrs,d.bcs(2)-d.bcs(1));
inds3 = ~(centers < 1 | centers > 39 | rs < 0.25 | PWs < 10);% | PWs > 20 | PWs < 10;

fitg = out4.fits.gauss; [rs,MFR,centers,PWs] = get_gauss_fit_parameters(fitg.coeffsrs,d.bcs(2)-d.bcs(1));
inds4 = ~(centers < 1 | centers > 39 | rs < 0.25 | PWs < 10);% | PWs > 20 | PWs < 10;

fitg = out5.fits.gauss; [rs,MFR,centers,PWs] = get_gauss_fit_parameters(fitg.coeffsrs,d.bcs(2)-d.bcs(1));
inds5 = ~(centers < 1 | centers > 39 | rs < 0.25 | PWs < 10);% | PWs > 20 | PWs < 10;
%%
zMIT_3 = props_C.zMI{4,1}; zMID_3 = props_C.zMI{4,7};
zMIT_4 = props_C.zMI{4,3}; zMID_4 = props_C.zMI{4,9};
%% visualize the data
out = out3;
while 1
    d.bcs = out.bin_centers;
    d.FR = out.FR_vs_speed;
    fitg = out.fits.gauss; fits = out.fits.sigmoid; fitl = out.fits.linear;
    d.fFRl = fitl.fitted; d.fFRs = fits.fitted; d.fFRg = fitg.fitted;
    d.cl = fitl.coeffsrs(:,3); d.cs = fits.coeffsrs(:,3); d.cg = fitg.coeffsrs(:,3);
    [rs,MFR,centers,PWs] = get_gauss_fit_parameters(fitg.coeffsrs,d.bcs(2)-d.bcs(1));
    inds = ~(centers < 1 | centers > 39 | rs < 0.25 | PWs < 10);% | PWs > 20 | PWs < 10;
    inds = inds3 | inds4 | inds5;
    % inds = logical(resp_speed{an,4});
%     t_resp = cell_list_op(resp,[],'or');
%     inds = ~inds;
%     inds = inds & ~t_resp{an}';
    100*sum(inds)/length(inds)
    d.FR = d.FR(inds,:); d.fFRl = d.fFRl(inds,:); d.fFRs = d.fFRs(inds,:); d.fFRg = d.fFRg(inds,:);
    d.cl = d.cl(inds); d.cs = d.cs(inds); d.cg = d.cg(inds);
    generic_multi_plot(1000,[3,4,size(d.FR,1)],'plotSpeedTuning',d)
    break;
end

%% confirmed here that the new uDataT variable with time binning has accurate values of firing rate as well as other variables
 % I compared the data in uDataT with the previously generated values of
 % Rs_C that I had with different types of events and then from there on
 % the different types of rasters etc.,
tRs = Rs_C{an,1};
ts = data_an.ts;
csel = data_an.air & data_an.C3;
figure(100);clf;plot(ts,csel);pause(0.05);
redges = find_rising_edge(csel,0.5,0.05);
fedges = find_falling_edge(csel,-0.5,0.05);
cell_num = 5;
cell_raster = []; lencr = [];
for ii = 1:length(redges)
    trial = redges(ii):fedges(ii);
    cell_raster{ii} = data_an.firing_rate(cell_num,trial);
    lencr(ii) = length(cell_raster{ii});
end
cell_raster1 = nan(10,max(lencr));
for ii = 1:length(redges)
    cell_raster1(ii,1:lencr(ii)) = cell_raster{ii};
end

trialnum = 10;
trial1_tRs = squeeze(tRs.sp_rasters(:,:,cell_num));
figure(100);clf;subplot(1,2,1);imagesc(trial1_tRs);colorbar;subplot(1,2,2);imagesc(cell_raster1);colorbar;

belt = find(data_an.belt);
data_an.ds(belt(4))-data_an.ds(belt(3))

%% now I was going to generate for each trial the variables ts, ds, sp, tm and then firing rate
% I combined the firing rates of the ON period (for example) and then
% concatenate them to see the overall mutual information.
% I will also find acceleration here
an = 1;
data_an = udataT{an};
% tRs = Rs_C{an,1};
ts = data_an.ts;ds = data_an.ds;sp = data_an.speed;tm = data_an.animal_motion;
ac = diff(sp)./diff(ts); ac = [0 ac];
csel = data_an.air & data_an.C5;
figure(100);clf;plot(ts,csel);pause(0.05);
redges = find_rising_edge(csel,0.5,0.05);
fedges = find_falling_edge(csel,-0.5,0.05);
cell_num = 14;
cell_raster = []; lencr = [];
atts = []; atds = []; atsp = []; attm = []; aFR = []; atac = [];
otts = [];
for ii = 1:length(redges)
    trial = redges(ii):fedges(ii);
    tts = ts(trial)-ts(trial(1)); tds = ds(trial)-ds(trial(1)); tsp = sp(trial); ttm = tm(trial); tac = ac(trial)
    FR = data_an.firing_rate(:,trial);
    atts = [atts tts]; atds = [atds tds]; atsp = [atsp tsp]; attm = [attm ttm]; aFR = [aFR FR]; atac = [atac tac];
    otts = [otts ts(trial)];
end
% find cells which have a lot of zeros
figure(100);clf;plot(FR');

firing_rate = data_an.firing_rate;

mFR = nanmean(firing_rate,2);
cells_to_keep = mFR > omFR(an);
sum(cells_to_keep)
% Create a new matrix excluding sparse cells
filtered_firing_rate = firing_rate(cells_to_keep, :);

%%
tic
nofmibins = 10; nshuffle = 500;
ctkaFR = aFR(cells_to_keep, :);
clear ts_O ds_O sp_O ac_O tm_O;
parfor jj = 1:size(ctkaFR,1)
    fr_neuron = ctkaFR(jj,:);
    [output ~] = info_metrics_S(fr_neuron, atts, nofmibins, [], nshuffle); ts_O(jj) = output;    
    [output ~] = info_metrics_S(fr_neuron, atds, nofmibins, [], nshuffle); ds_O(jj) = output;  
    [output ~] = info_metrics_S(fr_neuron, atsp, nofmibins, [], nshuffle); sp_O(jj) = output;
    [output ~] = info_metrics_S(fr_neuron, atac, nofmibins, [], nshuffle); ac_O(jj) = output;
    [output ~] = info_metrics_S(fr_neuron, attm, nofmibins, [], nshuffle); tm_O(jj) = output;
end
toc
beep
%%
ts_shp = arrayfun(@(s) s.('ShannonMI_p'), ts_O); ts_cp = arrayfun(@(s) s.('corr_p'), ts_O);

ds_shp = arrayfun(@(s) s.('ShannonMI_p'), ds_O); ds_cp = arrayfun(@(s) s.('corr_p'), ds_O);

sp_shp = arrayfun(@(s) s.('ShannonMI_p'), sp_O); sp_cp = arrayfun(@(s) s.('corr_p'), sp_O);

ac_shp = arrayfun(@(s) s.('ShannonMI_p'), ac_O); ac_cp = arrayfun(@(s) s.('corr_p'), ac_O);

tm_shp = arrayfun(@(s) s.('ShannonMI_p'), tm_O); tm_cp = arrayfun(@(s) s.('corr_p'), tm_O);

% sigvals = [
% [ts_shp < 0.05; ds_shp < 0.05; sp_shp < 0.05; ac_shp < 0.05; tm_shp < 0.05;]
% [ts_cp < 0.05; ds_cp < 0.05; sp_cp < 0.05; ac_cp < 0.05; tm_cp < 0.05;]
% ]

sigvals = [
[ts_shp < 0.05; ds_shp < 0.05; sp_shp < 0.05; ac_shp < 0.05;]
[ts_cp < 0.05; ds_cp < 0.05; sp_cp < 0.05; ac_cp < 0.05;]
]
ssigvals = sum(sigvals);
figure(100);clf;histogram(ssigvals, 'Normalization', 'probability')


% figure(100);clf;imagesc(sigvals)


%% visualize the graphs for the significant tuning
sel_var = 'ac';
cmdTxt = sprintf('selCells = ssigvals == 1 & (%s_cp < 0.05 | %s_shp < 0.05);',sel_var,sel_var)
eval(cmdTxt);
sum(selCells);
ctk = find(cells_to_keep);
idx_cells = find(selCells)
cn = 1;
%%
figure(100);clf;
% scatter(atts,ctkaFR(42,:))
plot(atts/max(atts));hold on;
plot(ctkaFR(43,:));
plot(atds/max(atds));
plot(atsp/max(atsp));
plot(atac/max(atac));


%%
nofmibins = 10; nshuffle = 0;
for jj = 1:size(FR,1)
    fr_neuron = aFR(jj,:);
    [output ~] = info_metrics_S(fr_neuron, atts, nofmibins, [], nshuffle);
    ts_MI(jj) = output.ShannonMI;
    [output ~] = info_metrics_S(fr_neuron, atds, nofmibins, [], nshuffle);
    ds_MI(jj) = output.ShannonMI;
    [output ~] = info_metrics_S(fr_neuron, atsp, nofmibins, [], nshuffle);
    sp_MI(jj) = output.ShannonMI;
    [output ~] = info_metrics_S(fr_neuron, attm, nofmibins, [], nshuffle);
    tm_MI(jj) = output.ShannonMI;
end
% [ts_MI(cell_num) ds_MI(cell_num) sp_MI(cell_num) tm_MI(cell_num)]
[ts_MI;ds_MI;sp_MI;tm_MI]
%%
figure(100);clf;plot(otts,aFR(14,:));