
clear all
clc
dataFolder = 'G:\Sam\OneDrive - University of Lethbridge\ExperimentalData\March 29, 2016, Tuesday'
fileName = '2016_03_29_0004.abf';
[d,si,h] = abf2load(makeName(fileName,dataFolder));

tiffDataFolder = 'G:\Sam\OneDrive - University of Lethbridge\ExperimentalData\March 29, 2016, Tuesday\Evoked_750FR_150Hz_Urethane_1mA_1pulse_LHP_0004';
tiffFileName = makeName('lo1.tif',tiffDataFolder);
[vsdi_image vsdi_image_info] = load_tiff_stack(tiffFileName);


% figure(1);clf;
% timeAxis = (1:h.dataPtsPerChan)*si*1e-6;
% plot(timeAxis,d(:,6,1));
% ylim([-0.001 3.5]);


% Get physical dimensions from first image because we know that our
% images' dimensions do not change in time.
% Our data was recorded with a frame rate of 150Hz, i.e. approx. one image every 6.7ms
% 5000 frames ~ 33.3s
FrameRate = 150;
width   = vsdi_image_info(1).Width;
height  = vsdi_image_info(1).Height;
nFrames = length(vsdi_image_info);
nPixels = width*height;

% Let's explore the data
ie_h = image_explorer(vsdi_image, vsdi_image_info, FrameRate, 'VSDI Spontaneous Activity - Unfiltered - RawData');


close all;

% Let's filter a bit and explore how things improve
% It makes live easier if we reshape our 3D image to 2D.
vsdi_vec = reshape(vsdi_image,nPixels,nFrames);

% Let's keep the memory clean. We do not need vsdi_image anymore.
% But let's keep the corresponding info.
clear vsdi_image;

% We later need the raw values of the first unfiltered frame to avoid some potential 
% numerical issues. Notice that it is just the first column in this
% representation.
vsdi_f0_offset = vsdi_vec(:,1);

% Let's bandpass filter !
%
% Recall: Aaron Gruber's great lecture on filters yesterday. We will use
% some of the items he talked about !

% There are many potential filters one can use, however, Majid's multi-year
% experience with many VSDI data sets has shown that a Type 1 Chebychev filter
% works well.
Nyquist     = FrameRate/2;
% You can play with these parameters; just rerun
% FROM HERE
FilterOrder = 2;             % was 2
RippleParam = 0.5;           % was 0.5
Bandpass    = [0.5, 6];      % was [0.5, 6]

% Construct a Chebychev Type 1 filter;
% Note: order n defines how many summands we have in the polynomials
% H(z) = (b(1)+b(2)z^{-1}+...+b(n-1)z^{-n}) /
% (1+a(2)z^{-1}+...+a(n-1)z^{-n}
[b,a] = cheby1(FilterOrder,RippleParam,Bandpass/Nyquist);

% If you like you can play with the filter parameters. Matlab provides a
% quick tool to study the filter's characteristics a.k.a. Bode plots.
figure;
freqz(b,a,Nyquist,FrameRate);

% Apply bandpass filter to our data
% WARNING:
% vsdi_vec_bp = filtfilt(b,a,double(vsdi_vec)); DOES NOT DO THE RIGHT JOB!
% filtfilt takes a vector and operates accross first non-singelton dim
vsdi_vec_filt = zeros(nPixels,nFrames);
for i=1:nPixels
    vsdi_vec_filt(i,:) = filtfilt(b,a,double(vsdi_vec(i,:)));
end


% Let's explore the results (Note: we can use the same info as before)
vsdi_image_filt = reshape(vsdi_vec_filt,width,height,nFrames);
ie_h = image_explorer(vsdi_image_filt, vsdi_image_info, FrameRate, 'VSDI Spontaneous Activity - Chebychev Type 1 Bandpass [0.5Hz-6Hz]');

% Chebychev Type 1 filters often introduce "edge artifacts" at the
% edges (beginning) which we want to remove. You can play with the cutoff
% if you like.
RemoveFront = 400;
RemoveBack  = 0;

% Valid Range
Start = RemoveFront+1;
End   = nFrames-RemoveBack;
vsdi_vec_filt_cut = vsdi_vec_filt(:,Start:End);
nFrames_cut = nFrames-(RemoveFront+RemoveBack);

vsdi_image_filt = reshape(vsdi_vec_filt_cut,width,height,nFrames_cut);
ie_h = image_explorer(vsdi_image_filt, vsdi_image_info(Start:End), FrameRate, 'VSDI Spontaneous Activity - Cutoff - Chebychev Type 1 Bandpass [0.5Hz-6Hz]');

% House keeping
vsdi_vec_filt = vsdi_vec_filt_cut;
nFrames = nFrames_cut;
vsdi_image_info = vsdi_image_info(Start:End);
clear Start End vsdi_vec_filt_cut nFrames_cut RemoveFront RemoveBack;

% Now let us normalize the data to obtain (f-f0)/f0 in percent
% After filtering our time signals are already quite flat and oscillate
% arround zero, i.e. their mean will be close to zero. Something close to
% zero is not a good choice for f0 because the quotiend diverges.
% Trick: Add the unfiltered values from the first frame to f before taking
% the mean to derive f0.
%ie_h = image_explorer(vsdi_image_filt, vsdi_image_info, FrameRate, 'VSDI Spontaneous Activity - Cutoff - Chebychev Type 1 Bandpass [0.5Hz-6Hz]');
vsdi_vec_filt = bsxfun(@plus,vsdi_vec_filt,vsdi_f0_offset);

% Do actual (f-f0)/f0 * 100% with f0 being the time-mean image
f0 = mean(vsdi_vec_filt,2);

for i=1:nFrames
    vsdi_vec_filt(:,i) = (vsdi_vec_filt(:,i)-f0)./f0 *100;
end

% And again, let's explore what we have
vsdi_image_filt = reshape(vsdi_vec_filt,width,height,nFrames);
ie_h = image_explorer(vsdi_image_filt, vsdi_image_info, FrameRate, 'VSDI Spontaneous Activity - Bandpass [0.5Hz-6Hz], df/f0');
