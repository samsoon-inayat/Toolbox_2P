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

% animal_list = {'NML_GC_01','NML_GC_01'};
% date_list = {'2025_12_15','2025_12_16'};
animal_list = {'NML_GC_01'};
date_list = {'2025_12_16'};

animal = get_exp_info(mD,animal_list,date_list);

disp('Done');
%%
colormaps = load('../../Common/Matlab/colorblind_colormap.mat');
colormaps.colorblind = flipud(colormaps.colorblind);
mD.colors = mat2cell(colormaps.colorblind,[ones(1,size(colormaps.colorblind,1))]);%{[0 0 0],[0.1 0.7 0.3],'r','b','m','c','g','y'}; % mD.colors = getColors(10,{'w','g'});
mD.dcolors = mat2cell(distinguishable_colors(20,'w'),[ones(1,20)]);
mD.axes_font_size = 6; mD.sigColor = [0.54 0.27 0.06];
mD.asterisk_font_size = 9;
mD.shades = generate_shades(3);
mD.pdf_folder = mD.adata_dir; 
mD.pd_folder = mD.pdata_dir;
mD.magfac = 1;
mData = mD;
disp('Done');

%%
if 1
%     make_db(T_C);
    owr = 0;
    animal = process_h264(animal,owr);
    animal = process_behavior_signals(animal);
end
disp('Done');

%%
% animal = air_signal_reader(animal);
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
plot(b.t,max(b.fSpeed)*(b.air_bin)/max(b.air_raw))
xlabel('Time (sec)');
ylabel('Speed cm/s');



%% Figure raw data
magfac = mD.magfac;
ff = makeFigureRowsCols(107,[3 5 6.5 1.5],'RowsCols',[1 1],'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.15 0.22],...
    'widthHeightAdjustment',[10 -200]);
MY = 2; ysp = 0.15285; mY = -2.5; titletxt = ''; ylabeltxt = {'PDF'}; % for all cells (vals) MY = 80
stp = 0.3*magfac; widths = [6.15 1 2.85 1]*magfac; gap = 0.115*magfac;
adjust_axes(ff,[mY MY],stp,widths,gap,{''});
axes_title_shifts_line = [0 0.55 0 0]; axes_title_shifts_text = [0.02 0.1 0 0]; xs_gaps = [1 2];
plot(b.tm,b.fSpeed);
hold on;
plot(b.tm,max(b.fSpeed)*(b.air_bin)/max(b.air_raw))
box off;
xlabel('Time (min)');
ylabel('Speed cm/s');
xlim([0 b.tm(end)]);
ylim([0 14.5])
format_axes(gca);

save_pdf(ff.hf,mD.pdf_folder,'bar_graph.pdf',600);  


%%
led_paws = extract_led_from_roi(b.led.paws, 60);

figure(100);clf;
plot(led_paws.time, led_paws.is_on); hold on
% stairs(led_paws.time, ...
%        double(led_paws.is_on) * max(led_paws.signal_s), ...
%        'LineWidth', 1.5)
xlabel('Time (s)')
ylabel('LED signal')
% legend('Smoothed ROI signal','LEfiD ON')
