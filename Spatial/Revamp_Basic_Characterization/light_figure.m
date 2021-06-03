function light_figure

mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size;
ei_11_15 = evalin('base','ei_11_15'); 
ei_2_3 = evalin('base','ei_2_3'); 

selContexts = [1 4 6];
rasterNames = {'light22T','light22T','light22T'};
Rs = get_rasters_data(ei_11_15,selContexts,rasterNames);
Rs = find_responsive_rasters(Rs,1:10);
n = 0;
%%
% traces and raster plots
for rr = 1:size(Rs,1)
    ccs = [];
    for cc = 1:size(Rs,2)
        R = Rs{rr,cc};
        ccs(:,cc) = R.resp.vals';
        resp_fraction(rr,cc) = R.resp.fraction;
    end
    resp_vals{rr} = ccs;
end
disp('done');
n = 0;

%%
an = 3; cn = 1;
plotRasters_simplest(Rs{an,cn})
