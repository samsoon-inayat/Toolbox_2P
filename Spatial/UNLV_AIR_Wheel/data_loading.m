for ii = 1:60
    pause(1);
end
%%
add_to_path_unlv
disp('Done')
%%
clear all
% clc
%%
main_dir = 'E:\GoogleDrive\InayatSamsoon\UNLV\AIR_Wheel_Methods'; mD.main_dir = main_dir;

rdata_dir = fullfile(main_dir,'RData'); mD.rdata_dir = rdata_dir;
pdata_dir = fullfile(main_dir,'PData'); mD.pdata_dir = pdata_dir;
adata_dir = fullfile(main_dir,'AData'); mD.adata_dir = adata_dir;

animal_list = {'NML_GC_01','NML_GC_01'};
date_list = {'2025_12_15','2025_12_16'};

animal = get_exp_info(mD,animal_list,date_list);

disp('Done');
%%
colormaps = load('../../../Common/Matlab/colorblind_colormap.mat');
colormaps.colorblind = flipud(colormaps.colorblind);
mD.colors = mat2cell(colormaps.colorblind,[ones(1,size(colormaps.colorblind,1))]);%{[0 0 0],[0.1 0.7 0.3],'r','b','m','c','g','y'}; % mD.colors = getColors(10,{'w','g'});
mD.dcolors = mat2cell(distinguishable_colors(20,'w'),[ones(1,20)]);
mD.axes_font_size = 6; mD.sigColor = [0.54 0.27 0.06];
mD.asterisk_font_size = 9;
mD.shades = generate_shades(3);
mD.pdf_folder = mD.adata_dir; 
mD.pd_folder = mD.pdata_dir;
mD.magfac = 1;
disp('Done');

%%
if 1
%     make_db(T_C);
    owr = 0;
    animal = process_h264(animal,owr);
    animal = process_behavior_signals(animal);
    for an = 1:length(animal)
        animal(an).b = load(animal(an).mat);
    end
end
disp('Done');

%%
animal = air_signal_reader(animal);
%%
b = animal(1).b;
b.dist = b.encoderCount * pi * 32/b.countsPerRev; % in cm
b.speed =  diff(b.dist)./diff(b.t); % in cm/sec
b.speed = double([0;b.speed]);
% b.speed = removeSpeedOutliers(b.speed);
% b.speed(b.speed < 0) = NaN;
% b.speed = fillmissing(b.speed,'linear',2,'EndValues','nearest');

samplingRate = 5000;
coeffs = ones(1, samplingRate)/samplingRate;
b.fSpeed = filter(coeffs, 1, b.speed);
%
figure(100);clf;
plot(b.t,b.fSpeed);
hold on;
plot(b.t,max(b.fSpeed)*(b.air_raw)/max(b.air_raw))