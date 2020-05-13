
%% Spectofram for Evoked Data

nFrames = 108;
FR = 150;
ImageSize = 128;

baseline_from = 20;
baseline_to = 30;

stimfrom = 1;
stimto = 5;

pathname = 'G:\Majid Comuter E\Data\Motif\2012\M072012_Awake\Evoked\VC_5ms_0.5%iso_post awake\';

file_hi = 'hi';
file_lo = 'lo';
file_no = 'no';
file_Mask = 'Mask.tif';

filename_hi = [pathname file_hi];
filename_lo = [pathname file_lo];
filename_no = [pathname file_no];
filename_Mask = [pathname file_Mask];

img_hi = single(zeros(ImageSize , ImageSize,nFrames));
img_lo = single(zeros(ImageSize , ImageSize,nFrames));
img_no = single(zeros(ImageSize , ImageSize,nFrames));
img_Mask = single(zeros(ImageSize , ImageSize));

img = [];
img_Mask  =  imreadalltiff(filename_Mask,1);
img_Mask = reshape(double(img_Mask),ImageSize * ImageSize,1);
[t1 t2] = find(img_Mask == 0);

for i = stimfrom:stimto
    
   filename1_hi = [filename_hi num2str(i) '.tif'] ;
   img_hi  =  single(imreadalltiff(filename1_hi,nFrames));
   
   filename1_lo = [filename_lo num2str(i) '.tif'] ;
   img_lo  =  single(imreadalltiff(filename1_lo,nFrames));
         
   filename1_no = [filename_no num2str(i) '.tif'] ;
   img_no  =  single(imreadalltiff(filename1_no,nFrames));
   
   img_hi = reshape (img_hi,ImageSize * ImageSize, nFrames);
   img_lo = reshape (img_lo,ImageSize * ImageSize, nFrames);
   img_no = reshape (img_no,ImageSize * ImageSize, nFrames);

      
   img_hi(:,31:33) = [];
   img_lo(:,31:33) = [];
   img_no(:,31:33) = [];
   
   
   img_hi = [img_hi(:,1:3) img_hi];
   img_lo = [img_lo(:,1:3) img_lo];
   img_no = [img_no(:,1:3) img_no];
   
   img_hi = reshape (img_hi,ImageSize , ImageSize, nFrames);
   img_lo = reshape (img_lo,ImageSize , ImageSize, nFrames);
   img_no = reshape (img_no,ImageSize , ImageSize, nFrames);
   
   [img_hi_out,img_lo_out]=trialbaseddff0(img_hi,img_lo,img_no,nFrames,baseline_from,baseline_to);

   img(1+(i-1)*2,:,:) = img_hi_out;
   img(2+(i-1)*2,:,:) = img_lo_out;
   
end

img(:,t1,:) = [];
img = reshape (img,size(img,1)*size(img,2),nFrames);

params.Fs = FR; % sampling frequency
fpassl=input('Give me the lower limit of the frequency range of interest between 0 and half of sampling frequency');
fpassu=input('Give me the upper limit of the frequency range of interest between 0 and half of sampling frequency. Upper limit must be > lower limit ');
if ~isempty(fpassl) && ~isempty(fpassu);
    if fpassl>fpassu; fpassl=0; display('Lower limit has to smaller than upper limit, using defaults, [0.1 50]'); end;
    if fpassu>params.Fs/2; fpassu=params.Fs/2;display('Upper limit has to smaller than half the sampling frequency,using defaults, [0.1 50]'); end;
else
    display('User did not specify the frequency range,using default, [0 50]');
    fpassl=0.1;
    fpassu=50;
end



defaultparams=0;
TW=input('Give me a time-bandwidth product >= 1 and < 10 ');
if ~isempty(TW);
    if TW>10; TW=1;display('Time-bandwidth product too high, using default, 1'); end;
    if TW<1; TW=1;display('Time-bandwidth product has to be greater than 1, using default, 1'); end;
else
    display('User has not specified time-bandwidth product, using default, 1');
    TW=1;
end

%TW = 2 or 3 yield the best results;
% mtspecgramc setting; 
params.tapers=[TW 2*TW-1];
params.fpass=[fpassl fpassu];
params.pad=6;
params.err=0;
params.trialave=1;


% movingwin=[0.075 0.0067] or movingwin=[0.1 0.0067] yields the best results;
twinl=input('Give me the moving window length for power spectrum between 0.01 and 2 in second ');
if isempty(twinl);
    display('User has not specified moving window length, using default, 0.5 s');
    twinl=0.5;
end
twinu=input('Give me the step window length for power spectrum between 0.01 and 0.5 second. moving window length must be > step window length ');
if isempty(twinu)
    display('User has not specified step window length, using default, 0.1 s');
    twinu=0.1;
end
movingwin=[twinl twinu];

clear S;
clear t;
clear f;

[S,t,f]=mtspecgramc(img',movingwin,params);


figure;
time = (0:1/FR:(nFrames/FR))';
time = time(1:end-1);
    
subplot(2,1,1),plot(time,mean(img),'Color',[0 0 1]);
set(gca,'FontName','Times New Roman','Fontsize', 16,'Position',[0.13 0.5838 0.7158 0.3132]);
ylabel('VSD amplitude (DF/F0(%))');
title('VSD Activity');
xlim([0 max(size(time))/FR])

position_figure
subplot(2,1,2),plot_matrix(S,t,f); xlabel([]); % plot spectrogram
caxis([-27.5 -24.75]); %colorbar; %Change the value to get the best contrast
set(gca,'FontName','Times New Roman','Fontsize', 16);
title({['VSD spectrogram,  W=' num2str(params.tapers(1)/movingwin(1)) 'Hz']; ['moving window = ' num2str(movingwin(1)) 's, step = ' num2str(movingwin(2)) 's']});
xlabel('time (s)'),ylabel('frequency Hz');
ylim([0 fpassu]);
xlim([0 max(size(time))/FR])

%%


X=getframe(gcf);
[file, pathname] = ...
    uigetfile('e:\Data\*.bmp','File Selector');
filenamelist=[pathname file];
file=filenamelist(1:end-4);
filefig=[file,'.bmp'];
if isempty(X.colormap)
    imwrite(X.cdata, filefig)
else
    imwrite(X.cdata, X.colormap, filefig)
end


