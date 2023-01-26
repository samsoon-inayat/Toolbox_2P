function o = get_data(ei,selContexts,rasterNames)

o.RsMC = get_rasters_dataMC(ei,selContexts,rasterNames); 
o.Rs = get_rasters_data(ei,selContexts,rasterNames); 
o.Rs = find_responsive_rasters(o.Rs,[]);
% o.Rs = find_responsive_rasters_1(o.Rs,[]);
% o.resp = get_responsive_cells(o.Rs); 
% [o.resp_FR,o.resp_FR_ei] = get_responsive_fraction_FR(ei,o.Rs,10000);
[o.mR,~,o.mR1] = calc_mean_rasters(o.Rs,[]); 
% o.props = get_props_Rs(o.Rs,50);

