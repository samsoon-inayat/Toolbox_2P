function plot_raster_efficient(fn,raster)

figure(fn);
xs = 1:size(raster,2);
nt = size(raster,1);
ys = reshape(raster',1,size(raster,2)*nt);

xs = repmat(xs,1,nt);

plot(xs,ys);


