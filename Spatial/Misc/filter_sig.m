function rasters = filter_sig(rasters)


n = 0;
Rs = rasters.sp_rasters1;



length_time_window = 1;
nsamples = ceil(length_time_window*FS); % secs
sigma = length_time_window/6;
w = gausswin(nsamples,(nsamples-1)/(2*sigma));
% for gausswin function sigma = (L-1)/(2a) from gausswin(L,a) ... so if filter time window is 0.3sec then sigma = 0.05 and if L = 3 then a = 20
%%
for ii = 1:size(spSigAll_i,1)
    tsig = spSigAll_i(ii,:);
    tsig = tsig + 1e-6;
    tsigf = filter(w,1,tsig);
    figure(1000);clf;plot(tsig);hold on;
    plot(tsigf);
    
end