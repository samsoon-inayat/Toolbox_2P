function fig_corr_map()


vp = evalin('base','vp');
vf = evalin('base','vf');
v = evalin('base','v');
mD = evalin('base','mData'); colors = mD.colors; sigColor = mD.sigColor; axes_font_size = mD.axes_font_size;
mData = mD;
animal = evalin('base','animal');

filen = fullfile(animal(1).pdir,'eye_pupil.mat');
if exist(filen,'file')
    load(filen)
end
n = 0;

%%
paws_energy_off = pre_means;
paws_energy_on = post_means;
%%
paws_OF_off = pre_means;
paws_OF_on = post_means;
%%


%%
face_energy_off = pre_means;
face_energy_on = post_means;
%%
face_OF_off = pre_means;
face_OF_on = post_means;
%%

%%
FR_speed_off = pre_means;
FR_speed_on = post_means;
%%
FL_speed_off = pre_means;
FL_speed_on = post_means;
%%
%%
HR_speed_off = pre_means;
HR_speed_on = post_means;
%%
HL_speed_off = pre_means;
HL_speed_on = post_means;
%%

%%
sface_energy_off = pre_means;
sface_energy_on = post_means;
%%
sface_OF_off = pre_means;
sface_OF_on = post_means;

%%
%%
eye_area_off = pre_means;
eye_area_on = post_means;

%%
%% =========================
%  Correlation maps: Air-ON vs Air-OFF (N=34 trials)
%  Assumes each variable below is a 34x1 vector.
%  =========================

ff = makeFigureRowsCols(2020,[5 4 6.9 3],'RowsCols',[1 3],'spaceRowsCols',[0 -0.001],'rightUpShifts',[0.081 0.2],'widthHeightAdjustment',[-20 -180]);


% ---- 1) Put variables in the order you want them displayed ----
varNames = { ...
    'PawsEnergy', 'PawsOF', ...
    'FaceEnergy', 'FaceOF', ...
    'FRspeed', 'FLspeed', 'HRspeed', 'HLspeed', ...
    'SideFaceEnergy', 'SideFaceOF', ...
    'EyeArea'};

% Build OFF and ON matrices (Ntrials x Nvars)
X_off = [ ...
    paws_energy_off(:), ...
    paws_OF_off(:), ...
    face_energy_off(:), ...
    face_OF_off(:), ...
    FR_speed_off(:), ...
    FL_speed_off(:), ...
    HR_speed_off(:), ...
    HL_speed_off(:), ...
    sface_energy_off(:), ...
    sface_OF_off(:), ...
    eye_area_off(:) ];

X_on  = [ ...
    paws_energy_on(:), ...
    paws_OF_on(:), ...
    face_energy_on(:), ...
    face_OF_on(:), ...
    FR_speed_on(:), ...
    FL_speed_on(:), ...
    HR_speed_on(:), ...
    HL_speed_on(:), ...
    sface_energy_on(:), ...
    sface_OF_on(:), ...
    eye_area_on(:) ];

% ---- 2) Sanity checks ----
nTrials_off = size(X_off,1);
nTrials_on  = size(X_on,1);
assert(nTrials_off == nTrials_on, 'ON and OFF must have same number of trials.');
assert(nTrials_on == 34, 'Expected 34 trials (but code will still run if not).');
assert(size(X_off,2) == numel(varNames), 'Mismatch: X_off columns vs varNames.');
assert(size(X_on,2)  == numel(varNames), 'Mismatch: X_on columns vs varNames.');

fprintf('Trials: %d | Vars: %d\n', nTrials_on, numel(varNames));
fprintf('NaNs (OFF): %d | NaNs (ON): %d\n', sum(isnan(X_off(:))), sum(isnan(X_on(:))));

% ---- 3) Correlation matrices (Pearson; pairwise to tolerate NaNs) ----
R_off = corr(X_off, 'Type','Pearson', 'Rows','pairwise'); R_off(eye(size(R_off)) == 1) = NaN;
R_on  = corr(X_on,  'Type','Pearson', 'Rows','pairwise'); R_on(eye(size(R_off)) == 1) = NaN;
dR    = R_on - R_off;  % ON minus OFF
dR(eye(size(R_off)) == 1) = NaN;

% ---- 4) Plot maps with consistent color limits ----
% For R matrices, standard [-1 1]. For dR, use symmetric limits.
dlim = max(abs(dR(:)));
if ~isfinite(dlim) || dlim == 0
    dlim = 1; % fallback
end

% figure('Color','w','Position',[100 100 1400 450]);

% (A) OFF
axes(ff.h_axes(1,1));% subplot(1,3,1);
imagesc(R_off, [-1 1]); axis image; colorbar;
title('Air-OFF');
set(gca,'XTick',1:numel(varNames),'XTickLabel',varNames,'XTickLabelRotation',45);
set(gca,'YTick',1:numel(varNames),'YTickLabel',varNames,'YTickLabelRotation',45);
% grid on;
box off;
set(gca,'Color','w')   % <-- THIS makes NaNs white
format_axes(gca)


% (B) ON
axes(ff.h_axes(1,2));%subplot(1,3,2);
imagesc(R_on, [-1 1]); axis image; colorbar;
title('Air-ON');
set(gca,'XTick',1:numel(varNames),'XTickLabel',varNames,'XTickLabelRotation',45);
set(gca,'YTick',1:numel(varNames),'YTickLabel',[]);
% grid on; 
box off;
set(gca,'Color','w')   % <-- THIS makes NaNs white
format_axes(gca)


% (C) Difference
axes(ff.h_axes(1,3));%subplot(1,3,3);
imagesc(dR, [-dlim dlim]); axis image; colorbar;
title('\DeltaR = ON - OFF');
set(gca,'XTick',1:numel(varNames),'XTickLabel',varNames,'XTickLabelRotation',45);
set(gca,'YTick',1:numel(varNames),'YTickLabel',[]);
box off;
set(gca,'Color','w')   % <-- THIS makes NaNs white
format_axes(gca)
save_pdf(ff.hf,mData.pdf_folder,sprintf('bar_graphs.pdf'),600);



% Optional: save
% exportgraphics(gcf, fullfile(mD.pdf_folder, 'Figure8_CorrMaps.pdf'), 'ContentType','vector');
