function raster_properties

mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size;
caSig_C = evalin('base','caSig_C'); caSig_A = evalin('base','caSig_A');
ei_C = evalin('base','ei10_C'); ei_A = evalin('base','ei10_A');
n = 0;
%%
bnd = [0.5 200];
for ani = 1:length(caSig_C)
    tcas = caSig_C{ani};
    tei = ei_C{ani};
    FR = tei.thorExp.frameRate;
    for cni = 1:size(tcas,1)
        cas = tcas(cni,:);
        [pL,fL,tL] = pspectrum(cas,FR,'spectrogram','FrequencyLimits',bnd,'TimeResolution',200);
        figure(2002);imagesc(tL,fL,pow2db(pL));set(gca,'Ydir','normal');colorbar
        pause(0.2);
    end
end
