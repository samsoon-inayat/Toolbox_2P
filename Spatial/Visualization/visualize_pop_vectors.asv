function visualize_pop_vectors

mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size;
ei = evalin('base','d15'); Rs = evalin('base','raster_data'); 
sC = evalin('base','selContexts'); rN = evalin('base','rasterNames');
an = 2;
cai_sampling_rate = ei{an}.thorExp.frameRate;
effective_sampling_rate = 1/0.15;
samplingRate = {'Ca','Ef','Ef','Ef','Ef','Ef','Ef','Ef','Ca','Ef'};
timeBefore = [2 5 7 NaN 7 NaN 7 NaN 2 5];

filename = fullfile(mData.pd_folder,'rasters.mat');

load(filename);

n = 0;
%%
ssC = [1 2 3 5 7 9 10];
for an = 1:5
    resp = [];
    for ii = 1:length(ssC)
        tR = rasters{an,ssC(ii)};
        resp(:,ii) = tR.resp.p < 0.05;
    end
    all_resp{an} = resp;
end
n = 0;
%%
sci = 3;
for an = 1:5
   resp = all_resp{an};
   tresp = resp(:,sci);
   tRs = get_selected_rasters(rasters{an,sci}.rasters,{tresp});
end
