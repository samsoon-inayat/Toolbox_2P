function MI = calc_MI(FR)
spatial_bin=1:size(FR,2); stimulus=repmat(spatial_bin,size(FR,1),1);
FRtemp=FR(:); FRtemp=FRtemp(~isnan(FRtemp));
stimulus=stimulus(:); stimulus=stimulus(~isnan(FR(:)));
FRbin=NaN(size(FRtemp)); FRbin1 = FRbin;
tempedges=[min(FRtemp) quantile(FRtemp(FRtemp>min(FRtemp)),linspace(0,1,4))];
tempedges(end)=tempedges(end)+0.01;
try
    [~,~,FRbin1] = histcounts(FRtemp,tempedges);
catch
end
MI = MutualInformation(stimulus,FRbin1);