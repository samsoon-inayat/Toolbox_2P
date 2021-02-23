function checking_ca_signal_from_within_roi
mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size;
ei = evalin('base','d15'); Rs = evalin('base','raster_data'); 
sC = evalin('base','selContexts'); rN = evalin('base','rasterNames');
 n = 0;
 %%
 figure(100);clf;
 an = 1; pp = 1;
 rec = ei{an};
 tP = rec.plane{pp};