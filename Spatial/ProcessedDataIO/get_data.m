function o = get_data(ei,selContexts,rasterNames)

o.Rs = get_rasters_data(ei,selContexts,rasterNames); 
o.Rs = find_responsive_rasters(o.Rs,[]);
o.resp = get_responsive_cells(o.Rs); 
o.resp_FR = get_responsive_fraction_FR(o.Rs);
o.mR = calc_mean_rasters(o.Rs,[]); 
o.props = get_props_Rs(o.Rs);

