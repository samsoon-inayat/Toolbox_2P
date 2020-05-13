clear all
clc

FrameRate = 150;
tiffDataFolder = 'G:\Sam\OneDrive - University of Lethbridge\ExperimentalData\March 29, 2016, Tuesday\Evoked_750FR_150Hz_Urethane_1mA_1pulse_LHP_0004';

fileName = sprintf('hi%d.tif',1);
tiffFileName = makeName(fileName,tiffDataFolder);
[vsdi_imageh vsdi_image_info] = load_tiff_stack(tiffFileName);


fileName = makeName('dFbyFo.mat',tiffDataFolder);
load(fileName);
ie_h = image_explorer(vsdi_image, vsdi_image_info, FrameRate, 'VSDI Activity');

tiffDataFolder = 'G:\Sam\OneDrive - University of Lethbridge\ExperimentalData\March 29, 2016, Tuesday\Evoked_750FR_150Hz_Urethane_1mA_100pulses_LHP_0005';
fileName = makeName('dFbyFo.mat',tiffDataFolder);
load(fileName);
ie_h = image_explorer(vsdi_image, vsdi_image_info, FrameRate, 'VSDI Activity');

tiffDataFolder = 'G:\Sam\OneDrive - University of Lethbridge\ExperimentalData\March 29, 2016, Tuesday\Evoked_750FR_150Hz_Urethane_1mA_1pulse_LHP_0006';
fileName = makeName('dFbyFo.mat',tiffDataFolder);
load(fileName);
ie_h = image_explorer(vsdi_image, vsdi_image_info, FrameRate, 'VSDI Activity');
