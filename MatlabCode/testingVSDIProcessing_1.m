%% clear variables in memory and clear console
clear all
clc


%% Load Data
% load ABF file
dataFolder = 'G:\Sam\OneDrive - University of Lethbridge\ExperimentalData\March 29, 2016, Tuesday'
fileName = '2016_03_29_0004.abf';
[d,si,h] = abf2load(makeName(fileName,dataFolder));

% load image stack
tiffDataFolder = 'G:\Sam\OneDrive - University of Lethbridge\ExperimentalData\March 29, 2016, Tuesday\Evoked_750FR_150Hz_Urethane_1mA_1pulse_LHP_0004';
tiffFileName = makeName('hi1.tif',tiffDataFolder);
[vsdi_image vsdi_image_info] = load_tiff_stack(tiffFileName);
FrameRate = 150;
width   = vsdi_image_info(1).Width;
height  = vsdi_image_info(1).Height;
nFrames = length(vsdi_image_info);
nPixels = width*height;
maskFileName = makeName('mask.tif','G:\Sam\OneDrive - University of Lethbridge\ExperimentalData\March 29, 2016, Tuesday\');
mask = imread(maskFileName);
for ii = 1:nFrames
   vsdi_image(:,:,ii) = vsdi_image(:,:,ii) .* mask;
end


%% view the 
ie_h = image_explorer(vsdi_image, vsdi_image_info, FrameRate, 'VSDI Spontaneous Activity - Unfiltered - RawData');

%% Reshape image data into fom 3D to 2D (vectorize)
vsdi_vec = reshape(vsdi_image,nPixels,nFrames);
vsdi_f0_offset = vsdi_vec(:,1);

%% Filter the data
Nyquist     = FrameRate/2;
FilterOrder = 2;             % was 2
RippleParam = 0.5;           % was 0.5
Bandpass    = [0.5, 6];      % was [0.5, 6]
[b,a] = cheby1(FilterOrder,RippleParam,Bandpass/Nyquist);
vsdi_vec_filt = zeros(nPixels,nFrames);
for i=1:nPixels
    vsdi_vec_filt(i,:) = filtfilt(b,a,double(vsdi_vec(i,:)));
end
% Let's explore the results (Note: we can use the same info as before)
vsdi_image_filt = reshape(vsdi_vec_filt,width,height,nFrames);
ie_h = image_explorer(vsdi_image_filt, vsdi_image_info, FrameRate, 'VSDI Spontaneous Activity - Chebychev Type 1 Bandpass [0.5Hz-6Hz]');


%% Chop off edges after chebysev filtering
RemoveFront = 50;
RemoveBack  = 0;
% Valid Range
Start = RemoveFront+1;
End   = nFrames-RemoveBack;
vsdi_vec_filt_cut = vsdi_vec_filt(:,Start:End);
nFrames_cut = nFrames-(RemoveFront+RemoveBack);
vsdi_image_filt = reshape(vsdi_vec_filt_cut,width,height,nFrames_cut);
ie_h = image_explorer(vsdi_image_filt, vsdi_image_info(Start:End), FrameRate, 'VSDI Spontaneous Activity - Cutoff - Chebychev Type 1 Bandpass [0.5Hz-6Hz]');

%% clearing some variables
vsdi_vec_filt = vsdi_vec_filt_cut;
nFrames = nFrames_cut;
vsdi_image_info = vsdi_image_info(Start:End);
clear Start End vsdi_vec_filt_cut nFrames_cut RemoveFront RemoveBack;


%% Find dF/Fo
% vsdi_vec_filt = vsdi_vec;
vsdi_vec_filt = bsxfun(@plus,vsdi_vec_filt,vsdi_f0_offset);
f0 = mean(vsdi_vec_filt,2);
for i=1:nFrames
    vsdi_vec_filt(:,i) = (vsdi_vec_filt(:,i)-f0)./f0 *100;
end
vsdi_image_filt = reshape(vsdi_vec_filt,width,height,nFrames);
ie_h = image_explorer(vsdi_image_filt, vsdi_image_info, FrameRate, 'VSDI Spontaneous Activity - Bandpass [0.5Hz-6Hz], df/f0');